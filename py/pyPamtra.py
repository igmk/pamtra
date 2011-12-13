# -*- coding: utf-8 -*-

#todo: function to calculate PIA?
from __future__ import division
import numpy as np
import datetime
import csv
import pickle
import time,calendar,datetime
import warnings
import sys
import os
from copy import deepcopy

import namelist #parser for fortran namelist files

try: 
	import numexpr as ne
except:
	warnings.warn("numexpr not available", Warning)

import meteoSI
try: 
	import pyPamtraLibWrapper
except:
	warnings.warn("pyPamtraLib not available", Warning)

try:
	import pp
except:
	warnings.warn("parallel python not available", Warning)

missingNumber=-9999

class pyPamtra(object):
	'''
	Class for pamtra calculations. Initalisations fill dictonary 'set' with default values and also fills 'dimensions' and 'units'
	
	'''
	def __init__(self):
		#set setting default values
		self.set = dict()
		self.set["verbose"]=0
		self.set["pyVerbose"] = 0 #not passed to Pamtra!
		self.set["dump_to_file"]=False
		self.set["tmp_path"]='/tmp/'
		self.set["data_path"]='/home/mech/models/pamtra/data/'

		self.set["obs_height"]=833000.
		self.set["units"]='T'
		self.set["outpol"]='VH'
		self.set["creator"]='pyPamtrauser'
		self.set["zeSplitUp"]=True # only locally in PYpamtra

		self.set["active"]=True
		self.set["passive"]=True

		self.set["ground_type"]='S'
		self.set["salinity"]=33.0
		self.set["EMissivity"]=0.6

		self.set["lgas_extinction"]=True
		self.set["gas_mod"]='R98'

		self.set["lhyd_extinction"]=True
		self.set["lphase_flag"]= True

		self.set["SD_cloud"]='C'
		
		self.set["SD_ice"]='C'
		self.set["EM_ice"]='mieic'
		
		self.set["SD_rain"]='C' 
		self.set["N_0rainD"]=8.0

		self.set["SD_snow"]='C' 
		self.set["N_0snowDsnow"]=7.628 
		self.set["EM_snow"]='densi' 
		self.set["snow_density"]=200 
		self.set["SP"]=0.2 
		self.set["isnow_n0"]=1
		self.set["liu_type"]=8

		self.set["SD_grau"]='C' 
		self.set["N_0grauDgrau"]=4.0 
		self.set["EM_grau"]='densi'
		self.set["graupel_density"]=400 

		self.set["SD_hail"]='C' 
		self.set["N_0hailDhail"]=4.0 
		self.set["EM_hail"]='densi'
		self.set["hail_density"]=917 

		self.set["n_moments"]=1
		self.set["moments_file"]='snowCRYSTAL'
		
		self.set["freqs"] = []
		self.set["nfreqs"] = 0
				
		self._setDefaultKeys = self.set.keys()
		
		self.dimensions = dict()
		
		self.dimensions["unixtime"] = ["ngridx","ngridy"]
		
		self.dimensions["ngridx"] = []
		self.dimensions["ngridy"] = []
		self.dimensions["nlyrs"] = ["ngridx","ngridy"]
		
		self.dimensions["model_i"] = ["ngridx","ngridy"]
		self.dimensions["model_j"] = ["ngridx","ngridy"]
		self.dimensions["lat"] = ["ngridx","ngridy"]
		self.dimensions["lon"] = ["ngridx","ngridy"]
		self.dimensions["lfrac"] = ["ngridx","ngridy"]
		self.dimensions["wind10u"] = ["ngridx","ngridy"]
		self.dimensions["wind10v"] = ["ngridx","ngridy"]
		
		self.dimensions["iwv"] = ["ngridx","ngridy"]
		self.dimensions["cwp"] = ["ngridx","ngridy"]
		self.dimensions["iwp"] = ["ngridx","ngridy"]
		self.dimensions["rwp"] = ["ngridx","ngridy"]
		self.dimensions["swp"] = ["ngridx","ngridy"]
		self.dimensions["gwp"] = ["ngridx","ngridy"]
		self.dimensions["hwp"] = ["ngridx","ngridy"]
		
		self.dimensions["hgt_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
		self.dimensions["temp_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
		self.dimensions["p_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
		self.dimensions["relhum_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
		
		self.dimensions["cwc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["iwc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["rwc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["swc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["gwc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["hwc_q"] = ["ngridx","ngridy","max_nlyrs"]
		
		self.dimensions["cwc_n"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["iwc_n"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["rwc_n"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["swc_n"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["gwc_n"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["hwc_n"] = ["ngridx","ngridy","max_nlyrs"]
		
		self.dimensions["hgt"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["Ze"] = ["gridx","gridy","lyr","frequency"]
		self.dimensions["Att_hydros"] = ["gridx","gridy","lyr","frequency"]
		self.dimensions["Att_atmo"] = ["gridx","gridy","lyr","frequency"]
		self.dimensions["tb"] = ["gridx","gridy","outlevels","angles","frequency","stokes"]

		
		self.units = dict()
		
		self.units["unixtime"] = "seconds since 1970-01-01 00:00"
		
		self.units["ngridx"] = "-"
		self.units["ngridy"] = "-"
		self.units["nlyrs"] = "-"
		
		self.units["model_i"] = "-"
		self.units["model_j"] = "-"
		self.units["lat"] = "deg.dec"
		self.units["lon"] = "deg.dec"
		self.units["lfrac"] = "-"
		self.units["wind10u"] = "m/s"
		self.units["wind10v"] = "m/s"
		
		self.units["iwv"] = "kg/m^2"
		self.units["cwp"] = "kg/m^2"
		self.units["iwp"] = "kg/m^2"
		self.units["rwp"] = "kg/m^2"
		self.units["swp"] = "kg/m^2"
		self.units["gwp"] = "kg/m^2"
		self.units["hwp"] = "kg/m^2"
		
		self.units["hgt_lev"] = "m"
		self.units["temp_lev"] = "K"
		self.units["p_lev"] = "Pa"
		self.units["relhum_lev"] = "1"
		
		self.units["cwc_q"] = "kg/kg"
		self.units["iwc_q"] = "kg/kg"
		self.units["rwc_q"] = "kg/kg"
		self.units["swc_q"] = "kg/kg"
		self.units["gwc_q"] = "kg/kg"
		self.units["hwc_q"] = "kg/kg"

		self.units["cwc_n"] = "#/kg"
		self.units["iwc_n"] = "#/kg"
		self.units["rwc_n"] = "#/kg"
		self.units["swc_n"] = "#/kg"
		self.units["gwc_n"] = "#/kg"
		self.units["hwc_n"] = "#/kg"
		
		self.units["hgt"] = "m"
		self.units["Ze"] = "dBz"
		self.units["Att_hydros"] = "-"
		self.units["Att_atmo"] = "-"
		self.units["tb"] = "K"
	
		self._nstokes = 2
		self._noutlevels = 2
		self._nangles = 32
		
		self.p = dict()
		self._helperP = dict()
		self.r = dict()
		
	def readNmlFile(self,inputFile):
		"""
		read classical Pamtra Namelist File from inputFile
		"""
		
		nmlFile = namelist.Namelist(inputFile)
		if nmlFile == {}:
			raise IOError("file not found: "+inputFile)
		
		for key in nmlFile.keys():
			for subkey in nmlFile[key]["par"][0].keys():
				if subkey in self._setDefaultKeys:
					if nmlFile[key]["par"][0][subkey][0] == ".true.": value = True
					elif nmlFile[key]["par"][0][subkey][0] == ".false.": value = False
					else: value = nmlFile[key]["par"][0][subkey][0]
					self.set[subkey] = value
					if self.set["pyVerbose"] > 1: print subkey, ":", value
				elif self.set["pyVerbose"] > 0:
					print "Setting '"+ subkey +"' from '"+ inputFile +"' skipped."
		if self.set["pyVerbose"] > 1: print "reading nml file done: ", inputFile
		return
		
	def readPamtraProfile(self,inputFile):
		"""
		read classical pamtra profile from file
		"""
		
		f = open(inputFile,"r")
		g = csv.reader(f,delimiter=" ",skipinitialspace=True)
		year,month,day,time, self.p["ngridx"], self.p["ngridy"], self.p["nlyrs"], self.p["deltax"], self.p["deltay"] = g.next()
		
		self.p["unixtime"] = calendar.timegm(datetime.datetime(year = int(year), month = int(month), day = int(day), hour = int(time[0:2]), minute = int(time[2:4]), second = 0).timetuple())
		
		self.p["ngridx"] = int(self.p["ngridx"])
		self.p["ngridy"] = int(self.p["ngridy"])
		self.p["nlyrs"] = int(self.p["nlyrs"])
		self.p["max_nlyrs"] = deepcopy(self.p["nlyrs"])
		
		if self.p["max_nlyrs"] > 200:
			warnings.warn("Too many layers for pamtra (max:200): " + str(self.p["max_nlyrs"]),Warning)
		
		
		self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
		self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],)
		self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"]+1,)
		
		self.p["model_i"] = np.zeros(self._shape2D)
		self.p["model_j"] = np.zeros(self._shape2D)
		self.p["lat"] = np.zeros(self._shape2D)
		self.p["lon"] = np.zeros(self._shape2D)
		self.p["lfrac"] = np.zeros(self._shape2D)
		self.p["wind10u"] = np.zeros(self._shape2D)
		self.p["wind10v"] = np.zeros(self._shape2D)
		
		self.p["iwv"] = np.zeros(self._shape2D)
		self.p["cwp"] = np.zeros(self._shape2D)
		self.p["iwp"] = np.zeros(self._shape2D)
		self.p["rwp"] = np.zeros(self._shape2D)
		self.p["swp"] = np.zeros(self._shape2D)
		self.p["gwp"] = np.zeros(self._shape2D)
		self.p["hwp"] = np.zeros(self._shape2D)
		
		self.p["hgt_lev"] = np.zeros(self._shape3Dplus)
		self.p["temp_lev"] = np.zeros(self._shape3Dplus)
		self.p["press_lev"] = np.zeros(self._shape3Dplus)
		self.p["relhum_lev"] = np.zeros(self._shape3Dplus)
		
		self.p["cwc_q"] = np.zeros(self._shape3D)
		self.p["iwc_q"] = np.zeros(self._shape3D)
		self.p["rwc_q"] = np.zeros(self._shape3D)
		self.p["swc_q"] = np.zeros(self._shape3D)
		self.p["gwc_q"] = np.zeros(self._shape3D)
		
		if self.set["n_moments"] == 1:
			self.p["hwc_q"] = np.ones(self._shape3D)*missingNumber
			self.p["cwc_n"] = np.ones(self._shape3D)*missingNumber
			self.p["iwc_n"] = np.ones(self._shape3D)*missingNumber
			self.p["rwc_n"] = np.ones(self._shape3D)*missingNumber
			self.p["swc_n"] = np.ones(self._shape3D)*missingNumber
			self.p["gwc_n"] = np.ones(self._shape3D)*missingNumber
			self.p["hwc_n"] = np.ones(self._shape3D)*missingNumber
			
			for x in np.arange(self.p["ngridx"]):
				for y in np.arange(self.p["ngridy"]):
					self.p["model_i"][x,y], self.p["model_j"][x,y] = np.array(np.array(g.next()),dtype=int)
					self.p["lat"][x,y], self.p["lon"][x,y], self.p["lfrac"][x,y],self.p["wind10u"][x,y],self.p["wind10v"][x,y]  = np.array(np.array(g.next()),dtype=float)
					self.p["iwv"][x,y],self.p["cwp"][x,y],self.p["iwp"][x,y],self.p["rwp"][x,y],self.p["swp"][x,y],self.p["gwp"][x,y] = np.array(np.array(g.next()),dtype=float)
					self.p["hgt_lev"][x,y,0],self.p["press_lev"][x,y,0],self.p["temp_lev"][x,y,0],self.p["relhum_lev"][x,y,0] = np.array(np.array(g.next()),dtype=float)
					for z in np.arange(self.p["nlyrs"]):
						self.p["hgt_lev"][x,y,z+1],self.p["press_lev"][x,y,z+1],self.p["temp_lev"][x,y,z+1],self.p["relhum_lev"][x,y,z+1],self.p["cwc_q"][x,y,z],self.p["iwc_q"][x,y,z],self.p["rwc_q"][x,y,z],self.p["swc_q"][x,y,z],self.p["gwc_q"][x,y,z] = np.array(np.array(g.next()),dtype=float)

		elif self.set["n_moments"] == 2:
			self.p["hwc_q"] = np.zeros(self._shape3D)
			self.p["cwc_n"] = np.zeros(self._shape3D)
			self.p["iwc_n"] = np.zeros(self._shape3D)
			self.p["rwc_n"] = np.zeros(self._shape3D)
			self.p["swc_n"] = np.zeros(self._shape3D)
			self.p["gwc_n"] = np.zeros(self._shape3D)
			self.p["hwc_n"] = np.zeros(self._shape3D)
			
			for x in np.arange(self.p["ngridx"]):
				for y in np.arange(self.p["ngridy"]):
					self.p["model_i"][x,y], self.p["model_j"][x,y] = np.array(np.array(g.next()),dtype=int)
					self.p["lat"][x,y], self.p["lon"][x,y], self.p["lfrac"][x,y],self.p["wind10u"][x,y],self.p["wind10v"][x,y]  = np.array(np.array(g.next()),dtype=float)
					self.p["iwv"][x,y],self.p["cwp"][x,y],self.p["iwp"][x,y],self.p["rwp"][x,y],self.p["swp"][x,y],self.p["gwp"][x,y],self.p["hwp"][x,y] = np.array(np.array(g.next()),dtype=float)
					self.p["hgt_lev"][x,y,0],self.p["press_lev"][x,y,0],self.p["temp_lev"][x,y,0],self.p["relhum_lev"][x,y,0] = np.array(np.array(g.next()),dtype=float)
					for z in np.arange(self.p["nlyrs"]):
						self.p["hgt_lev"][x,y,z+1],self.p["press_lev"][x,y,z+1],self.p["temp_lev"][x,y,z+1],self.p["relhum_lev"][x,y,z+1],self.p["cwc_q"][x,y,z],self.p["iwc_q"][x,y,z],self.p["rwc_q"][x,y,z],self.p["swc_q"][x,y,z],self.p["gwc_q"][x,y,z],self.p["hwc_q"][x,y,z],self.p["cwc_n"][x,y,z],self.p["iwc_n"][x,y,z],self.p["rwc_n"][x,y,z],self.p["swc_n"][x,y,z],self.p["gwc_n"][x,y,z],self.p["hwc_n"][x,y,z] = np.array(np.array(g.next()),dtype=float)
		else:
			raise IOError("n_moments mus be 1 or 2")


		#in PyPamtra I want relhum not in %
		self.p["relhum_lev"] = self.p["relhum_lev"]/100.

		#finally make an array from nlyrs and unixtime
		self.p["nlyrs"] = np.ones(self._shape2D,dtype=int)*self.p["nlyrs"]
		self.p["unixtime"] = np.ones(self._shape2D,dtype=int)*self.p["unixtime"]

		f.close()

	def createProfile(self,**kwargs):
		'''
		Function to create the Pamtra Profiles.
		
		Variables ending on _lev mean level values, variables without are layer values (height vecotr one entry shorter!)!
		
		Everything is needed in SI units, relhum is in Pa/PA not %
		
		The following variables are mandatroy:
		hgt_lev, temp_lev, press_lev and (relhum_lev OR q)
		
		The following variables are optional and guessed if not provided:	"timestamp","lat","lon","lfrac","wind10u","wind10v","hgt_lev","cwc_q","iwc_q","rwc_q","swc_q","gwc_q"
		
		The integrated values are calculated if not provided:
		"iwv","cwp","iwp","rwp","swp","gwp"
		
		The following variables are only needed together with the 2 moments scheme!
		"hwc_q","hwp", "cwc_n","iwc_n","rwc_n","swc_n","gwc_n","hwc_n"
		'''
		allVars=["timestamp","lat","lon","lfrac","wind10u","wind10v","iwv","cwp","iwp","rwp","swp","gwp","hwp","hgt_lev","press_lev","temp_lev","relhum_lev","q","cwc_q","iwc_q","rwc_q","swc_q","gwc_q","hwc_q","cwc_n","iwc_n","rwc_n","swc_n","gwc_n","hwc_n"]
			
		for key in kwargs.keys():
			if key not in allVars:
				raise TypeError("Could not parse "+key)
		
		if not ("hgt_lev" in kwargs.keys() and "temp_lev" in kwargs.keys() and "press_lev" in kwargs.keys() and ("relhum_lev" in kwargs.keys() or "q" in kwargs.keys())):
			raise TypeError("I need hgt_lev and temp_lev and press_lev and (relhum_lev or q)!")
		
		noDims = len(np.shape(kwargs["hgt_lev"]))
		
		if noDims == 1:
			self.p["ngridx"] = 1
			self.p["ngridy"] = 1
		elif noDims == 2:
			self.p["ngridx"] = np.shape(kwargs["hgt_lev"])[0]
			self.p["ngridy"] = 1
		elif noDims == 3:
			self.p["ngridx"] = np.shape(kwargs["hgt_lev"])[0]
			self.p["ngridy"] = np.shape(kwargs["hgt_lev"])[1]
		else:
			print "Too many dimensions!"
			raise IOError
		
		self.p["max_nlyrs"] = np.shape(kwargs["hgt_lev"])[-1] -1
		
		if self.p["max_nlyrs"] > 200:
			warnings.warn("Too many layers for pamtra (max:200): " + str(self.p["max_nlyrs"]),Warning)
		
		#if np.any(self.p["nlyrs"] != self.p["max_nlyrs"]):
			#self._radiosonde = True
		#else:
			#self._radiosonde = False
			
		self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
		self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
		self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
		
		self.p["nlyrs"] = np.sum(kwargs["hgt_lev"]!=missingNumber,axis=-1) -1
		self.p["nlyrs"] = self.p["nlyrs"].reshape(self._shape2D)
		
		self.p["hgt_lev"] = kwargs["hgt_lev"].reshape(self._shape3Dplus)
		self.p["temp_lev"] = kwargs["temp_lev"].reshape(self._shape3Dplus)
		self.p["press_lev"] = kwargs["press_lev"].reshape(self._shape3Dplus)

		self.p["deltax"] = 0.
		self.p["deltay"] = 0.

		self.p["model_i"] = np.array(np.where(self.p["press_lev"][:,:,0])[0]).reshape(self._shape2D) +1
		self.p["model_j"] = np.array(np.where(self.p["press_lev"][:,:,0])[1]).reshape(self._shape2D) +1

		if "timestamp" not in kwargs.keys():
			self.p["unixtime"] = np.ones(self._shape2D)* int(time.time())
			warnings.warn("timestamp set to now", Warning)
		else:
			if type(kwargs["timestamp"]) in (int,np.int32,np.int64,float,np.float32,np.float64) :
				self.p["unixtime"] = np.ones(self._shape2D,dtype=int)*kwargs["timestamp"]
			elif (type(kwargs["timestamp"]) == np.ndarray):
				if kwargs["timestamp"].dtype in (int,np.int32,np.int64,float,np.float32,np.float64):
					self.p["unixtime"] = kwargs["timestamp"].reshape(self._shape2D)
				else:
					raise TypeError("timestamp entries have to be int or float objects")
			elif (type(kwargs["timestamp"]) == datetime):
				self.p["unixtime"] = np.ones(self._shape2D,dtype=int)*calendar.timegm(kwargs["timestamp"].timetuple())
			else:
				raise TypeError("timestamp has to be int, float or datetime object")
			
		for environment, preset in [["lat",50.938056],["lon",6.956944],["lfrac",1],["wind10u",0],["wind10v",0]]:
			if environment not in kwargs.keys():
				self.p[environment] = np.ones(self._shape2D)*preset
				warnings.warn("%s set to %s"%(environment,preset,), Warning)
			else:
				if type(kwargs[environment]) in (int,np.int32,np.int64,float,np.float32,np.float64):
					self.p[environment] = np.ones(self._shape2D) * kwargs[environment]
				else:
					self.p[environment] = kwargs[environment].reshape(self._shape2D)
		
		for qValue in ["cwc_q","iwc_q","rwc_q","swc_q","gwc_q","hwc_q"]:
			if qValue not in kwargs.keys():
				self.p[qValue] = np.zeros(self._shape3D)
				warnings.warn(qValue + " set to 0", Warning)
			else:
				self.p[qValue] = kwargs[qValue].reshape(self._shape3D)
				
		for hydroNo in ["cwc_n","iwc_n","rwc_n","swc_n","gwc_n","hwc_n"]:
			if hydroNo not in kwargs.keys():
				self.p[hydroNo] = np.ones(self._shape3D) * -9999
				warnings.warn(hydroNo + " set to -9999", Warning)
			else:
				self.p[hydroNo] = kwargs[hydroNo].reshape(self._shape3D)


		if "q" in kwargs.keys():
			self._helperP["q"] = kwargs["q"].reshape(self._shape3D)

		if "relhum_lev" in kwargs.keys():
			self.p["relhum_lev"] = kwargs["relhum_lev"].reshape(self._shape3Dplus)
		else:
			self._calcRelhum_lev()

		for pDict,qValue,intValue in [[self._helperP,"q","iwv"],[self.p,"cwc_q","cwp"],[self.p,"iwc_q","iwp"],[self.p,"rwc_q","rwp"],[self.p,"swc_q","swp"],[self.p,"gwc_q","gwp"],[self.p,"hwc_q","hwp"]]:
			if intValue in kwargs.keys():
				self.p[intValue] = kwargs[intValue].reshape(self._shape2D)
			else:
				#now we need q!
				self._calcQ()
				#nothing to do without hydrometeors:
				if np.all(pDict[qValue]==0):
					self.p[intValue] = np.zeros(self._shape2D)
				else:
					self._calcMoistRho() #provides also temp,press,dz and q!
					self.p[intValue] =  np.sum(pDict[qValue]*self._helperP["rho_moist"]*self._helperP["dz"],axis=-1)
					
		#clean up: remove all nans
		for key in self.p.keys():
			#apply only to arrays, lists or tupels
			if np.sum(np.shape(self.p[key])) > 0:
				self.p[key][np.isnan(self.p[key])] = missingNumber
		
		return

	def _calcPressTempRelhum(self):
		if set(("temp","relhum","press","dz")).issubset(self._helperP.keys()):
			return
		else:
			if self.set["pyVerbose"] > 0:
				print 'calculating "temp","relhum","press","dz"'
			self._helperP["temp"] = (self.p["temp_lev"][...,0:-1] + self.p["temp_lev"][...,1:])/2.
			self._helperP["relhum"] = (self.p["relhum_lev"][...,0:-1] + self.p["relhum_lev"][...,1:])/2.
			dz = np.diff(self.p["hgt_lev"],axis=-1)
			self._helperP["dz"] = dz
			
			#if not self._radiosonde:

				#p0 = self.p["press_lev"][...,0:-1]
				#p1 = self.p["press_lev"][...,1:]
				#xp = ne.evaluate("-1.*log(p1/p0)/dz")

			#else: #to allow array operations, cases with varying heights or invalid data have to be treated a bit differently:

			dz[dz<=0]=9999
			self._helperP["relhum"][self._helperP["relhum"]<=0] = 1
			self._helperP["temp"][self._helperP["temp"]<=0] = 1
			
			press_lev1 = deepcopy(self.p["press_lev"])
			press_lev1[press_lev1==missingNumber]=1
			
			p0 = press_lev1[...,0:-1]
			p1 = press_lev1[...,1:]
			xp = ne.evaluate("-1.*log(p1/p0)/dz")
			xp[xp==0] = 9999
				
			self._helperP["press"] = ne.evaluate("-1.*p0/xp*(exp(-xp*dz)-1.)/dz")
			del p0,p1,xp
				
			return

	def _calcQ(self):
		
		if set(("q",)).issubset(self._helperP.keys()):
			return
		else:
			self._calcPressTempRelhum()
			
			if self.set["pyVerbose"] > 0:
				print 'calculating "q"'
			self._helperP["q"] = meteoSI.rh2q(self._helperP["relhum"],self._helperP["temp"],self._helperP["press"])
			return

	def _calcQ_lev(self):
		if set(("q_lev",)).issubset(self._helperP.keys()):
			return
		else:
			if self.set["pyVerbose"] > 0:
				print 'calculating "q_lev"'
			qBot = self._helperP["q"][...,0:1] + 0.25*(self._helperP["q"][...,0:1]-self._helperP["q"][...,1:2])
			qTop = self._helperP["q"][...,-1:] + 0.25*(self._helperP["q"][...,-1:]-self._helperP["q"][...,-2:-1])
			qMid = (self._helperP["q"][...,0:-1] + self._helperP["q"][...,1:])/2.
			
			self._helperP["q_lev"] = np.concatenate((qBot,qMid,qTop),axis=-1)
			return

	def _calcMoistRho(self):
		if set(("rho_moist",)).issubset(self._helperP.keys()):
			return
		else:
			self._calcQ()
			self._calcPressTempRelhum()
			
			if self.set["pyVerbose"] > 0:
				print 'calculating "rho_moist"'

			self._helperP["rho_moist"] = meteoSI.moist_rho_q(self._helperP["press"],self._helperP["temp"],self._helperP["q"],self.p["cwc_q"],self.p["iwc_q"],self.p["rwc_q"],self.p["gwc_q"],self.p["swc_q"])
			return

	def _calcRelhum_lev(self):
		if set(("relhum_lev",)).issubset(self.p.keys()):
			return
		else:
			self._calcQ_lev()
			
			if self.set["pyVerbose"] > 0:
				print 'calculating "relhum_lev"'
				
			self.p["relhum_lev"] = meteoSI.q2rh(self._helperP["q_lev"],self.p["temp_lev"],self.p["press_lev"])
			return




	def filterProfiles(self,condition):
		'''
		discard profiles, which do not fullfill "condition" (2D boolean array)
		
		Note: If the initial field is a 2D grid, the field is flattend after application
		
		Applicate between CreateProfile and runPamtra
		
		'''
		if condition.shape != self._shape2D and condition.shape != (self._shape2D[0]*self._shape2D[1],):
			raise ValueError("shape mismatch, shape of condition must be 2D field!")
		
		condition = condition.reshape(self._shape2D)
		
		#create a new shape!
		self.p["ngridx"] = np.sum(condition)
		self.p["ngridy"] = 1
		
		self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
		self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
		self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)

		for key in ["unixtime","nlyrs","lat","lon","lfrac","model_i","model_j","wind10u","wind10v","iwv","cwp","iwp","rwp","swp","gwp","hwp"]:
			self.p[key] = self.p[key][condition].reshape(self._shape2D)

		for key in ["cwc_q","iwc_q","rwc_q","swc_q","gwc_q","hwc_q","cwc_n","iwc_n","rwc_n","swc_n","gwc_n","hwc_n"]:
			self.p[key] = self.p[key][condition].reshape(self._shape3D)
			
		for key in ["hgt_lev","temp_lev","press_lev","relhum_lev"]:
			self.p[key] = self.p[key][condition].reshape(self._shape3Dplus)
		
		for key in ['q_lev', 'temp', 'q', 'dz', 'press', 'relhum', 'rho_moist']:
			if key in self._helperP.keys():
				if key in ['q_lev']:
					self._helperP[key] = self._helperP[key][condition].reshape(self._shape3Dplus)
				if key in ['temp', 'q', 'dz', 'press', 'relhum', 'rho_moist']:
					self._helperP[key] = self._helperP[key][condition].reshape(self._shape3D)
		return

	def rescaleHeights(self,new_hgt_lev):
	
		# sort height vectors
		old_hgt_lev = self.p["hgt_lev"]
		new_hgt = (new_hgt_lev[...,1:] + new_hgt_lev[...,:-1])/2.
		old_hgt = (old_hgt_lev[...,1:] + old_hgt_lev[...,:-1])/2.
		
		self.p["max_nlyrs"] = len(new_hgt_lev) -1
		
		#new shape
		self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
		self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
		

		for key in ["cwc_q","iwc_q","rwc_q","swc_q","gwc_q","hwc_q","cwc_n","iwc_n","rwc_n","swc_n","gwc_n","hwc_n"]:
			#make new array
			newP = np.ones(self._shape3D) * missingNumber
			for x in xrange(self._shape2D[0]):
				for y in xrange(self._shape2D[1]):
					#interpolate!
					newP[x,y] = np.interp(new_hgt,old_hgt[x,y],self.p[key][x,y])
			#save new array
			self.p[key] = newP
			#and mark all entries below -1 as missing Number!
			self.p[key][self.p[key]<-1] = missingNumber
			
		for key in ["hgt_lev","temp_lev","relhum_lev"]:
			newP = np.ones(self._shape3Dplus) * missingNumber
			for x in xrange(self._shape2D[0]):
				for y in xrange(self._shape2D[1]):
					newP[x,y] = np.interp(new_hgt_lev,old_hgt_lev[x,y],self.p[key][x,y])
			self.p[key] = newP
			self.p[key][self.p[key]<-1] = missingNumber
			
		for key in ["press_lev"]:
			newP = np.ones(self._shape3Dplus) * missingNumber
			for x in xrange(self._shape2D[0]):
				for y in xrange(self._shape2D[1]):
					newP[x,y] = np.exp(np.interp(new_hgt_lev,old_hgt_lev[x,y],np.log(self.p[key][x,y])))
			self.p[key] = newP
			self.p[key][self.p[key]<-1] = missingNumber

		self.p["nlyrs"] = np.sum(self.p["hgt_lev"] != missingNumber,axis=-1) -1
		
		if self.p["max_nlyrs"] > 200:
			warnings.warn("Still too many layers for pamtra (max:200): " + str(self.p["max_nlyrs"]),Warning)


		return

	def addCloudShape(self):
		"""
		adds cloud base and cloud top to the data to an existing pamtra profile
		clouds are detected with rh >= 0.95 and T >= 253.15 K
		if multiple cloud layers are found, only the uppermost and the lowest one are saved.
		"""
		self.p["cloudTop"] = np.ones(self._shape2D)*missingNumber
		self.p["cloudBase"] = np.ones(self._shape2D)*missingNumber
		
		for x in xrange(self._shape2D[0]):
			for y in xrange(self._shape2D[1]):
				i_top, i_base, i_cloud =  meteoSI.detect_liq_cloud(self.p["hgt_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["temp_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["relhum_lev"][x,y,0:self.p["nlyrs"][x,y]])
				if len(i_top)> 0:
					self.p["cloudTop"][x,y] = self.p["hgt_lev"][x,y,i_top[-1]+1]
					self.p["cloudBase"][x,y] = self.p["hgt_lev"][x,y,i_base[0]]
		return

	def addPseudoAdiabaticLWC(self):
		"""
		adds liquid water content to the data to an existing pamtra profile
		using a simple cloud model with pseudo adiabatic lapse rate
		clouds are detected with rh >= 0.95 and T >= 253.15 K
		existing cwc_q values are overriden!
		"""

		#calculate CWC
		self.p["cwc_q"][:] = 0

		for x in xrange(self._shape2D[0]):
			for y in xrange(self._shape2D[1]):
				i_top, i_base, i_cloud =  meteoSI.detect_liq_cloud(self.p["hgt_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["temp_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["relhum_lev"][x,y,0:self.p["nlyrs"][x,y]])
				for k in np.arange(len(i_top)):
					i_cloud_k = np.arange(i_base[k],i_top[k]+1) 
					self.p["cwc_q"][x,y,i_cloud_k[:-1]] = meteoSI.mod_ad(self.p["temp_lev"][x,y,i_cloud_k], self.p["press_lev"][x,y,i_cloud_k], self.p["hgt_lev"][x,y,i_cloud_k], 1)

		self.p["cwp"][:] = 0
		
		self._calcMoistRho() #provides also temp,press,dz and q!
		self.p["cwp"][:] =  np.sum(self.p["cwc_q"]*self._helperP["rho_moist"]*self._helperP["dz"],axis=-1)

		return

	def createFullProfile(self,timestamp,lat,lon,lfrac,wind10u,wind10v,
			iwv,cwp,iwp,rwp,swp,gwp,hwp,
			hgt_lev,press_lev,temp_lev,relhum_lev,
			cwc_q,iwc_q,rwc_q,swc_q,gwc_q,
			hwc_q,cwc_n,iwc_n,rwc_n,swc_n,gwc_n,hwc_n):
		
		'''
		create comple PAmtra Profile
		
		No Extras, no missing values are guessed. the Data is only reshaped
		'''
		
		
		noDims = len(np.shape(hgt_lev))
		
		if noDims == 1:
			self.p["ngridx"] = 1
			self.p["ngridy"] = 1
		elif noDims == 2:
			self.p["ngridx"] = np.shape(hgt_lev)[0]
			self.p["ngridy"] = 1
		elif noDims == 3:
			self.p["ngridx"] = np.shape(hgt_lev)[0]
			self.p["ngridy"] = np.shape(hgt_lev)[1]
		else:
			print "Too many dimensions!"
			raise IOError
		
		self.p["nlyrs"] = np.sum(hgt_lev!=missingNumber,axis=-1) -1
		self.p["max_nlyrs"] = np.shape(hgt_lev)[-1] -1
		
		if self.p["max_nlyrs"] > 200:
			warnings.warn("Too many layers for pamtra (max:200): " + str(self.p["max_nlyrs"]),Warning)
		
		self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
		self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
		self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
		

		self.p["unixtime"] = timestamp.reshape(self._shape2D)
				
		self.p["deltax"] = 0.
		self.p["deltay"] = 0.
		self.p["lat"] = lat.reshape(self._shape2D)
		self.p["lon"] = lon.reshape(self._shape2D)
		self.p["lfrac"] = lfrac.reshape(self._shape2D)
		self.p["model_i"] = np.array(np.where(lat.reshape(self._shape2D))[0]).reshape(self._shape2D) +1
		self.p["model_j"] = np.array(np.where(lon.reshape(self._shape2D))[1]).reshape(self._shape2D) +1
		self.p["wind10u"] = wind10u.reshape(self._shape2D)
		self.p["wind10v"] = wind10v.reshape(self._shape2D)

		self.p["iwv"] = iwv.reshape(self._shape2D)
		self.p["cwp"] = cwp.reshape(self._shape2D)
		self.p["iwp"] = iwp.reshape(self._shape2D)
		self.p["rwp"] = rwp.reshape(self._shape2D)
		self.p["swp"] = swp.reshape(self._shape2D)
		self.p["gwp"] = gwp.reshape(self._shape2D)
		self.p["hwp"] = gwp.reshape(self._shape2D)
		
		self.p["hgt_lev"] = hgt_lev.reshape(self._shape3Dplus)
		self.p["temp_lev"] = temp_lev.reshape(self._shape3Dplus)
		self.p["press_lev"] = press_lev.reshape(self._shape3Dplus)
		self.p["relhum_lev"] = relhum_lev.reshape(self._shape3Dplus)
		
		self.p["cwc_q"] = cwc_q.reshape(self._shape3D)
		self.p["iwc_q"] = iwc_q.reshape(self._shape3D)
		self.p["rwc_q"] = rwc_q.reshape(self._shape3D)
		self.p["swc_q"] = swc_q.reshape(self._shape3D)
		self.p["gwc_q"] = gwc_q.reshape(self._shape3D)
		
		self.p["hwc_q"] = hwc_q.reshape(self._shape3D)
		self.p["cwc_n"] = cwc_n.reshape(self._shape3D)
		self.p["iwc_n"] = iwc_n.reshape(self._shape3D)
		self.p["rwc_n"] = rwc_n.reshape(self._shape3D)
		self.p["swc_n"] = swc_n.reshape(self._shape3D)
		self.p["gwc_n"] = gwc_n.reshape(self._shape3D)
		self.p["hwc_n"] = hwc_n.reshape(self._shape3D)

	def runParallelPamtra(self,freqs,pp_servers=(),pp_local_workers="auto",pp_deltaF=0,pp_deltaX=0,pp_deltaY = 0):
		'''
		run Pamtra analouge to runPamtra, but with with parallel python
		
		input:
		
		freqs: list with frequencies
		pp_servers: tuple with servers, ("*") activates auto detection
		pp_local_workers: number of local workers, "auto" is usually amount of cores
		pp_deltaF,pp_deltaX,pp_deltaY: length of the peaces in frequency,x and y direction. 0 means no slicing
		
		In my experience, smaller pieces (e.g. pp_deltaF=0,pp_deltaX=1,pp_deltaY=1) work much better than bigger ones, even though overhead might be bigger.
		
		output is collected by _ppCallback()
		'''
		
		tttt = time.time()
		
		if np.max(self.p["relhum_lev"])>5:
			raise IOError("relative humidity is _not_ in %!")

		
		for key in self.set:
			if key not in self._setDefaultKeys:
				warnings.warn("Warning can not parse setting: ",key, Warning)

		if pp_local_workers == "auto":
			self.job_server = pp.Server(ppservers=pp_servers,secret="pyPamtra") 
		else:
			self.job_server = pp.Server(pp_local_workers,ppservers=pp_servers,secret="pyPamtra") 
			
			
		if self.set["pyVerbose"] >= 0: 
			print "Starting pp with: "
			pp_nodes = self.job_server.get_active_nodes()
			for key in pp_nodes.keys():
				print key+": "+str(pp_nodes[key])+" nodes"
		
		if type(freqs) in (int,np.int32,np.int64,float,np.float32,np.float64): freqs = [freqs]
		
		self.set["freqs"] = freqs
		self.set["nfreqs"] = len(freqs)
		
		if pp_deltaF==0: pp_deltaF = self.set["nfreqs"]
		if pp_deltaX==0: pp_deltaX = self.p["ngridx"]
		if pp_deltaY==0: pp_deltaY = self.p["ngridy"]
		
		pp_ii = -1
		pp_jobs = dict()
		
		self.r["Ze"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Att_hydro"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Att_atmo"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		
		self.r["Ze_cw"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Ze_rr"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Ze_ci"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Ze_sn"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Ze_gr"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Ze_ha"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Att_cw"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Att_rr"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Att_ci"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Att_sn"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Att_gr"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
		self.r["Att_ha"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber

		self.r["hgt"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],))*missingNumber
		self.r["tb"] = np.ones((self.p["ngridx"],self.p["ngridy"],self._noutlevels,self._nangles,self.set["nfreqs"],self._nstokes))*missingNumber
		
		self.pp_noJobs = len(np.arange(0,self.set["nfreqs"],pp_deltaF))*len(np.arange(0,self.p["ngridx"],pp_deltaX))*len(np.arange(0,self.p["ngridy"],pp_deltaY))
		self.pp_jobsDone = 0
		
		fi = open("/tmp/pp_logfile.txt","w")
		fi.write("Starting pp with %i jobs \n\r"%(self.pp_noJobs))
		fi.close()
		
		self.hosts=[]
		
		for pp_startF in np.arange(0,self.set["nfreqs"],pp_deltaF):
			pp_endF = pp_startF + pp_deltaF
			if pp_endF > self.set["nfreqs"]: pp_endF = self.set["nfreqs"]
			pp_nfreqs = pp_endF - pp_startF
			for pp_startX in np.arange(0,self.p["ngridx"],pp_deltaX):
				pp_endX = pp_startX + pp_deltaX
				if pp_endX > self.p["ngridx"]: pp_endX = self.p["ngridx"]
				pp_ngridx = pp_endX - pp_startX
				for pp_startY in np.arange(0,self.p["ngridy"],pp_deltaY):
					pp_endY = pp_startY + pp_deltaY
					if pp_endY > self.p["ngridy"]: pp_endY = self.p["ngridy"]
					pp_ngridy = pp_endY - pp_startY
					
					pp_ii+=1
					
					pp_jobs[pp_ii] = self.job_server.submit(pyPamtraLibWrapper.PamtraFortranWrapper, (
					#self.set
					self.set["verbose"], self.set["dump_to_file"], self.set["tmp_path"], self.set["data_path"], self.set["obs_height"], self.set["units"], self.set["outpol"], self.set["creator"], self.set["active"], self.set["passive"], self.set["ground_type"], self.set["salinity"], self.set["EMissivity"], self.set["lgas_extinction"], self.set["gas_mod"], self.set["lhyd_extinction"], self.set["lphase_flag"],self.set["SD_cloud"], self.set["SD_ice"], self.set["EM_ice"], self.set["SD_rain"], self.set["N_0rainD"], self.set["SD_snow"], self.set["N_0snowDsnow"], self.set["EM_snow"], self.set["snow_density"], self.set["SP"], self.set["isnow_n0"], self.set["liu_type"], self.set["SD_grau"], self.set["N_0grauDgrau"], self.set["EM_grau"], self.set["graupel_density"], self.set["SD_hail"], self.set["N_0hailDhail"], self.set["EM_hail"], self.set["hail_density"], self.set["n_moments"], self.set["moments_file"],
					#input
					pp_ngridx,
					pp_ngridy,
					self.p["max_nlyrs"],
					self.p["nlyrs"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					pp_nfreqs,
					self.set["freqs"][pp_startF:pp_endF],
					self.p["unixtime"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["deltax"],self.p["deltay"],
					self.p["lat"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["lon"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["model_i"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["model_j"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["wind10u"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["wind10v"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["lfrac"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["relhum_lev"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["press_lev"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["temp_lev"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["hgt_lev"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["iwv"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["cwp"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["iwp"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["rwp"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["swp"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["gwp"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["hwp"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["cwc_q"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["iwc_q"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["rwc_q"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["swc_q"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["gwc_q"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["hwc_q"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["cwc_n"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["iwc_n"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["rwc_n"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["swc_n"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["gwc_n"][pp_startX:pp_endX,pp_startY:pp_endY].tolist(),
					self.p["hwc_n"][pp_startX:pp_endX,pp_startY:pp_endY].tolist()
					),tuple(), ("pyPamtraLibWrapper","pyPamtraLib","os",),callback=self._ppCallback,
					callbackargs=(pp_startX,pp_endX,pp_startY,pp_endY,pp_startF,pp_endF,pp_ii,))
					
					if self.set["pyVerbose"] >= 0: 
						sys.stdout.write("\r"+20*" "+"\r"+ "%i, %5.3f%% submitted"%(pp_ii+1,(pp_ii+1)/float(self.pp_noJobs)*100))
						sys.stdout.flush()

		if self.set["pyVerbose"] >= 0: 
			print " "
			print self.pp_noJobs, "jobs submitted"


		if self.set["pyVerbose"] >= 0: print " "; self.job_server.get_active_nodes()
		self.job_server.wait()
		
		self.r["settings"] = self.set
		self.r["pamtraVersion"] = self.r["pamtraVersion"].strip()
		self.r["pamtraHash"] = self.r["pamtraHash"].strip()
		
		if self.set["pyVerbose"] >= 0: print " "; self.job_server.print_stats()
		self.job_server.destroy()
		del self.job_server
		
		if self.set["pyVerbose"] >= 0: print "pyPamtra runtime:", time.time() - tttt
	
	def _ppCallback(self,pp_startX,pp_endX,pp_startY,pp_endY,pp_startF,pp_endF,pp_ii,*results):
		'''
		Collect the data of parallel pyPamtra		
		'''
		(((
		self.r["pamtraVersion"],self.r["pamtraHash"], 
		self.r["Ze"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Ze_cw"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Ze_rr"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Ze_ci"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Ze_sn"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Ze_gr"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Ze_ha"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Att_hydro"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Att_cw"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Att_rr"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Att_ci"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Att_sn"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Att_gr"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Att_ha"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["Att_atmo"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF], 
		self.r["hgt"][pp_startX:pp_endX,pp_startY:pp_endY,:],
		self.r["tb"][pp_startX:pp_endX,pp_startY:pp_endY,:,:,pp_startF:pp_endF,:], 
		self.r["angles"],
		),host,),) = results
		self.pp_jobsDone += 1
		if self.set["pyVerbose"] > 0: 
			sys.stdout.write("\r"+50*" "+"\r"+ "%s: %6i, %8.3f%% collected (#%6i, %s)"%(datetime.datetime.now().strftime("%Y%m%d-%H:%M:%S"),self.pp_jobsDone,(self.pp_jobsDone)/float(self.pp_noJobs)*100,pp_ii+1,host))
			sys.stdout.flush()
		if self.set["pyVerbose"] > 1: print " "; self.job_server.print_stats()
		fi = open("/tmp/pp_logfile.txt","a")
		fi.write("%s: %6i, %8.3f%% collected (#%6i, %s)\n\r"%(datetime.datetime.now().strftime("%Y%m%d-%H:%M:%S"),self.pp_jobsDone,(self.pp_jobsDone)/float(self.pp_noJobs)*100,pp_ii+1,host))
		fi.close()
		
	def runPamtra(self,freqs):
		'''
		run Pamtra from python
		'''
		tttt = time.time()
		
		if type(freqs) == int in (int,np.int32,np.int64,float,np.float32,np.float64): freqs = [freqs]
		
		self.set["freqs"] = freqs
		self.set["nfreqs"] = len(freqs)
		
		for key in self.set:
			if key not in self._setDefaultKeys:
				print "Warning Could not parse ",key
		
		if np.max(self.p["relhum_lev"])>5:
			raise IOError("relative humidity is _not_ in %!")

		#output
		(
		self.r["pamtraVersion"],
		self.r["pamtraHash"],
		self.r["Ze"], 
		self.r["Ze_cw"], 
		self.r["Ze_rr"], 
		self.r["Ze_ci"], 
		self.r["Ze_sn"], 
		self.r["Ze_gr"], 
		self.r["Ze_ha"], 
		self.r["Att_hydro"], 
		self.r["Att_cw"], 
		self.r["Att_rr"], 
		self.r["Att_ci"], 
		self.r["Att_sn"], 
		self.r["Att_gr"], 
		self.r["Att_ha"], 
		self.r["Att_atmo"], 
		self.r["hgt"],
		self.r["tb"],
		self.r["angles"],
		), host = \
		pyPamtraLibWrapper.PamtraFortranWrapper(
		#self.set
		self.set["verbose"], self.set["dump_to_file"], self.set["tmp_path"], self.set["data_path"], self.set["obs_height"], self.set["units"], self.set["outpol"], self.set["creator"], self.set["active"], self.set["passive"], self.set["ground_type"], self.set["salinity"], self.set["EMissivity"], self.set["lgas_extinction"], self.set["gas_mod"], self.set["lhyd_extinction"], self.set["lphase_flag"],self.set["SD_cloud"], self.set["SD_ice"], self.set["EM_ice"], self.set["SD_rain"], self.set["N_0rainD"], self.set["SD_snow"], self.set["N_0snowDsnow"], self.set["EM_snow"], self.set["snow_density"], self.set["SP"], self.set["isnow_n0"], self.set["liu_type"], self.set["SD_grau"], self.set["N_0grauDgrau"], self.set["EM_grau"], self.set["graupel_density"], self.set["SD_hail"], self.set["N_0hailDhail"], self.set["EM_hail"], self.set["hail_density"], self.set["n_moments"], self.set["moments_file"],
		#input
		self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.p["nlyrs"],self.set["nfreqs"],self.set["freqs"],
		self.p["unixtime"],
		self.p["deltax"],self.p["deltay"], self.p["lat"],self.p["lon"],self.p["model_i"],self.p["model_j"],
		self.p["wind10u"],self.p["wind10v"],self.p["lfrac"],
		self.p["relhum_lev"],self.p["press_lev"],self.p["temp_lev"],self.p["hgt_lev"],
		self.p["iwv"],self.p["cwp"],self.p["iwp"],self.p["rwp"],self.p["swp"],self.p["gwp"],self.p["hwp"],
		self.p["cwc_q"],self.p["iwc_q"],self.p["rwc_q"],self.p["swc_q"],self.p["gwc_q"],
		self.p["hwc_q"],self.p["cwc_n"],self.p["iwc_n"],self.p["rwc_n"],self.p["swc_n"],self.p["gwc_n"],self.p["hwc_n"])
				
		self.r["settings"] = self.set
		
		self.r["pamtraVersion"] = self.r["pamtraVersion"].strip()
		self.r["pamtraHash"] = self.r["pamtraHash"].strip()
		
		#for key in self.__dict__.keys():
			#print key
			#print self.__dict__[key]
		
		if self.set["pyVerbose"] >= 0: print "pyPamtra runtime:", time.time() - tttt


	def writeResultsToNumpy(self,fname):
		'''
		write the complete state of the session (profile,results,settings to a file
		'''
		f = open(fname, "w")
		pickle.dump([self.r,self.p,self.set], f)
		f.close()

	def loadResultsFromNumpy(self,fname):
		'''
		load the complete state of the session (profile,results,settings to a file
		'''
		try: 
			f = open(fname, "r")
			[self.r,self.p,self.set] = pickle.load(f)
			f.close()
		except:
			raise IOError ("Could not read data")
		
	def writeResultsToNetCDF(self,fname,profileVars="all",ncForm="NETCDF3_CLASSIC"):
		'''
		write the results to a netcdf file
		
		Input:
		
		fname: str filename with path
		profileVars list of variables of the profile to be saved. "all" saves all implmented ones
		ncForm: str netcdf file format, possible values are NETCDF3_CLASSIC, NETCDF3_64BIT, NETCDF4_CLASSIC, and NETCDF4
		'''
		import netCDF4
		try: 
			self.r
			self.r["pamtraVersion"]
		except:
			raise IOError ("run runPamtra first!")

		cdfFile = netCDF4.Dataset(fname,"w",ncForm)
		
		#write meta data
		cdfFile.history = "Created with pyPamtra (Version: "+self.r["pamtraVersion"]+", Git Hash: "+self.r["pamtraHash"]+") by "+self.set["creator"]+" (University of Cologne, IGMK) at " + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
		
		cdfFile.properties = str(self.set)

		#make frequsnions
		cdfFile.createDimension('grid_x',self.p["ngridx"])
		cdfFile.createDimension('grid_y',self.p["ngridy"])
		cdfFile.createDimension('frequency',self.set["nfreqs"])
		
		if (self.r["settings"]["passive"]):
			cdfFile.createDimension('angles',len(self.r["angles"]))
			cdfFile.createDimension('outlevels',self._noutlevels)
			cdfFile.createDimension('stokes',self._nstokes)
			
		if (self.r["settings"]["active"]):
			cdfFile.createDimension('heightbins',self.p["max_nlyrs"])
		
		dim2d = ("grid_x","grid_y",)
		dim3d = ("grid_x","grid_y","heightbins",)
		dim4d = ("grid_x","grid_y","heightbins","frequency")
		dim6d = ("grid_x","grid_y","outlevels","angles","frequency","stokes")
			
		#create and write dim variables
		
		nc_frequency = cdfFile.createVariable('frequency','f4',('frequency',),fill_value= missingNumber)
		nc_frequency.units = 'GHz'
		nc_frequency[:] = self.set["freqs"]
		
		nc_gridx = cdfFile.createVariable('grid_x','f4',('grid_x',),fill_value= missingNumber)
		nc_gridx.units = '-'
		nc_gridx[:] = np.arange(self.p["ngridx"])

		nc_gridy = cdfFile.createVariable('grid_y','f4',('grid_y',),fill_value= missingNumber)
		nc_gridy.units = '-'
		nc_gridy[:] = np.arange(self.p["ngridy"])

		
		
		if (self.r["settings"]["active"]):

			nc_heightbins = cdfFile.createVariable('heightbins', 'i4',("heightbins",),fill_value= missingNumber)
			nc_heightbins.units = "-"
			nc_heightbins[:] = np.arange(0,self.p["max_nlyrs"])


			nc_height = cdfFile.createVariable('height', 'f4',dim3d,fill_value= missingNumber)
			nc_height.units = "m"
			nc_height[:] = self.r["hgt"]
		
		if (self.r["settings"]["passive"]):
			nc_angle = cdfFile.createVariable('angles','f4',('angles',),fill_value= missingNumber)
			nc_angle.units = 'deg'
			nc_angle[:] = self.r["angles"]
			
			nc_stokes = cdfFile.createVariable('stokes', 'S1',("stokes",),fill_value= missingNumber)
			nc_stokes.units = "-"
			nc_stokes[:] = ["H","V"]
			
			nc_out = cdfFile.createVariable('outlevels', 'f4',("outlevels",),fill_value= missingNumber)
			nc_out.units = "m over sea level (top of atmosphere value) OR m over ground (ground value)"
			nc_out[:] = [self.set["obs_height"],0]
			
		#create and write variables
		
		nc_model_i = cdfFile.createVariable('model_i', 'i4',dim2d,fill_value= missingNumber)
		nc_model_i.units = "-"
		nc_model_i[:] = self.p["model_i"]
		
		nc_model_j = cdfFile.createVariable('model_j', 'i4',dim2d,fill_value= missingNumber)
		nc_model_j.units = "-"
		nc_model_j[:] = self.p["model_j"]
		
		nc_nlyrs = cdfFile.createVariable('nlyr', 'i4',dim2d,fill_value= missingNumber)
		nc_nlyrs.units = "-"
		nc_nlyrs[:] = self.p["nlyrs"]
		
		nc_time = cdfFile.createVariable('datatime', 'i4',dim2d,fill_value= missingNumber)
		nc_time.units = "seconds since 1970-01-01 00:00:00"
		nc_time[:] = self.p["unixtime"]
		
		nc_longitude = cdfFile.createVariable('longitude', 'f4',dim2d,fill_value= missingNumber)
		nc_longitude.units = "deg.dec"
		nc_longitude[:] = self.p["lon"]
		
		nc_latitude = cdfFile.createVariable('latitude', 'f4',dim2d,fill_value= missingNumber)
		nc_latitude.units = "deg.dec"
		nc_latitude[:] = self.p["lat"]
		
		nc_lfrac = cdfFile.createVariable('lfrac', 'f4',dim2d,fill_value= missingNumber)
		nc_lfrac.units = "-"
		nc_lfrac[:] = self.p["lfrac"]
		
		
		if (self.r["settings"]["active"]):
			
			if self.set["zeSplitUp"]:

				nc_Ze_cw = cdfFile.createVariable('Ze_cloud_water', 'f4',dim4d,fill_value= missingNumber)
				nc_Ze_cw.units = "dBz"
				nc_Ze_cw[:] = self.r["Ze_cw"]

				nc_Ze_rr = cdfFile.createVariable('Ze_rain', 'f4',dim4d,fill_value= missingNumber)
				nc_Ze_rr.units = "dBz"
				nc_Ze_rr[:] = self.r["Ze_rr"]

				nc_Ze_ci = cdfFile.createVariable('Ze_cloud_ice', 'f4',dim4d,fill_value= missingNumber)
				nc_Ze_ci.units = "dBz"
				nc_Ze_ci[:] = self.r["Ze_ci"]

				nc_Ze_sn = cdfFile.createVariable('Ze_snow', 'f4',dim4d,fill_value= missingNumber)
				nc_Ze_sn.units = "dBz"
				nc_Ze_sn[:] = self.r["Ze_sn"]

				nc_Ze_gr = cdfFile.createVariable('Ze_graupel', 'f4',dim4d,fill_value= missingNumber)
				nc_Ze_gr.units = "dBz"
				nc_Ze_gr[:] = self.r["Ze_gr"]

				nc_Att_cw = cdfFile.createVariable('Attenuation_cloud_water', 'f4',dim4d,fill_value= missingNumber)
				nc_Att_cw.units = "dB"
				nc_Att_cw[:] = self.r["Att_cw"]

				nc_Att_rrrs = cdfFile.createVariable('Attenuation_rain', 'f4',dim4d,fill_value= missingNumber)
				nc_Att_rrrs.units = "dB"
				nc_Att_rrrs[:] = self.r["Att_rr"]
				
				nc_Att_ci = cdfFile.createVariable('Attenuation_cloud_ice', 'f4',dim4d,fill_value= missingNumber)
				nc_Att_ci.units = "dB"
				nc_Att_ci[:] = self.r["Att_ci"]
				
				nc_Att_sn = cdfFile.createVariable('Attenuation_snow', 'f4',dim4d,fill_value= missingNumber)
				nc_Att_sn.units = "dB"
				nc_Att_sn[:] = self.r["Att_sn"]
				
				nc_Att_gr = cdfFile.createVariable('Attenuation_graupel', 'f4',dim4d,fill_value= missingNumber)
				nc_Att_gr.units = "dB"
				nc_Att_gr[:] = self.r["Att_gr"]
				
				if self.set["n_moments"]==2:

					nc_Ze_ha = cdfFile.createVariable('Ze_hail', 'f4',dim4d,fill_value= missingNumber)
					nc_Ze_ha.units = "dBz"
					nc_Ze_ha[:] = self.r["Ze_ha"]

					nc_Att_ha = cdfFile.createVariable('Attenuation_hail', 'f4',dim4d,fill_value= missingNumber)
					nc_Att_ha.units = "dB"
					nc_Att_ha[:] = self.r["Att_ha"]
					
				
			else: 
				 
				nc_Ze = cdfFile.createVariable('Ze', 'f4',dim4d,fill_value= missingNumber)
				nc_Ze.units = "dBz"
				nc_Ze[:] = self.r["Ze"]
				
				nc_Attenuation_Hydrometeors = cdfFile.createVariable('Attenuation_Hydrometeors', 'f4',dim4d,fill_value= missingNumber)
				nc_Attenuation_Hydrometeors.units = "dB"
				nc_Attenuation_Hydrometeors[:] = self.r["Att_hydro"]
			
			nc_Attenuation_Atmosphere = cdfFile.createVariable('Attenuation_Atmosphere', 'f4',dim4d,fill_value= missingNumber)
			nc_Attenuation_Atmosphere.units = "dB"
			nc_Attenuation_Atmosphere[:] = self.r["Att_atmo"]
		
		if (self.r["settings"]["passive"]):
			nc_tb = cdfFile.createVariable('tb', 'f4',dim6d,fill_value= missingNumber)
			nc_tb.units = "K"
			nc_tb[:] = self.r["tb"]
			
		
		#profile data
		
		if "iwv" in profileVars or profileVars =="all":
			nc_iwv = cdfFile.createVariable('iwv', 'f4',dim2d,fill_value= missingNumber)
			nc_iwv.units = "kg/m^2"
			nc_iwv[:] = self.p["iwv"]

		if "cwp" in profileVars or profileVars =="all":
			nc_cwp= cdfFile.createVariable('cwp', 'f4',dim2d,fill_value= missingNumber)
			nc_cwp.units = "kg/m^2"
			nc_cwp[:] = self.p["cwp"]
			
		if "iwp" in profileVars or profileVars =="all":
			nc_iwp = cdfFile.createVariable('iwp', 'f4',dim2d,fill_value= missingNumber)
			nc_iwp.units = "kg/m^2"
			nc_iwp[:] = self.p["iwp"]

		if "rwp" in profileVars or profileVars =="all":
			nc_rwp = cdfFile.createVariable('rwp', 'f4',dim2d,fill_value= missingNumber)
			nc_rwp.units = "kg/m^2"
			nc_rwp[:] = self.p["rwp"]
			
		if "swp" in profileVars or profileVars =="all":
			nc_swp = cdfFile.createVariable('swp', 'f4',dim2d,fill_value= missingNumber)
			nc_swp.units = "kg/m^2"
			nc_swp[:] = self.p["swp"]
			
		if "gwp" in profileVars or profileVars =="all":
			nc_gwp = cdfFile.createVariable('gwp', 'f4',dim2d,fill_value= missingNumber)
			nc_gwp.units = "kg/m^2"
			nc_gwp[:] = self.p["gwp"]

		if "hwp" in profileVars or profileVars =="all":
			nc_hwp = cdfFile.createVariable('hwp', 'f4',dim2d,fill_value= missingNumber)
			nc_hwp.units = "kg/m^2"
			nc_hwp[:] = self.p["hwp"]
			
		if ("cloudBase" in profileVars or profileVars =="all") and ("cloudBase" in self.p.keys()):
			nc_cb = cdfFile.createVariable('cloudBase', 'f4',dim2d,fill_value= missingNumber)
			nc_cb.units = "m"
			nc_cb[:] = self.p["cloudBase"]
			
		if ("cloudTop" in profileVars or profileVars =="all") and ("cloudTop" in self.p.keys()):
			nc_ct = cdfFile.createVariable('cloudTop', 'f4',dim2d,fill_value= missingNumber)
			nc_ct.units = "m"
			nc_ct[:] = self.p["cloudTop"]

		cdfFile.close()
		if self.set["pyVerbose"] >= 0: print fname,"written"



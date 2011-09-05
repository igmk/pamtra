# -*- coding: utf-8 -*-

import numpy as np
import datetime
from pyPamtraLib import pypamtralib
import pyPamtraLib
import csv
import pickle
import time,calendar, datetime

def PamtraFortranWrapper(set_verbose, set_dump_to_file, set_tmp_path, set_data_path, set_obs_height, set_units, set_outpol, set_creator, set_active, set_passive, set_ground_type, set_salinity, set_emissivity, set_lgas_extinction, set_gas_mod, set_lhyd_extinction, set_lphase_flag, set_SD_snow, set_N_0snowDsnow, set_EM_snow, set_SP, set_isnow_n0, set_liu_type, set_SD_grau, set_N_0grauDgrau, set_EM_grau, set_EM_ice, set_SD_rain, set_N_0rainD, set_n_moments, set_moments_file,
	#input
	profile_nlyr,profile_ngridx,profile_ngridy,nfreqs,freqs,
	profile_year,profile_month,profile_day,profile_time,
	profile_deltax,profile_deltay, profile_lat,profile_lon,profile_model_i,profile_model_j,
	profile_wind10u,profile_wind10v,profile_lfrac,
	profile_relhum_lev,profile_press_lev,profile_temp_lev,profile_hgt_lev,
	profile_iwv,profile_cwp,profile_iwp,profile_rwp,profile_swp,profile_gwp,
	profile_cwc_q,profile_iwc_q,profile_rwc_q,profile_swc_q,profile_gwc_q):
	
	
	result_pamtraVersion, result_pamtraHash,\
	result_Ze, result_attenuationHydro, result_attenuationAtmo, result_hgt, result_tb, result_angles = \
	pypamtralib(
	#self.set
	set_verbose, set_dump_to_file, set_tmp_path, set_data_path, set_obs_height, set_units, set_outpol, set_creator, set_active, set_passive, set_ground_type, set_salinity, set_emissivity, set_lgas_extinction, set_gas_mod, set_lhyd_extinction, set_lphase_flag, set_SD_snow, set_N_0snowDsnow, set_EM_snow, set_SP, set_isnow_n0, set_liu_type, set_SD_grau, set_N_0grauDgrau, set_EM_grau, set_EM_ice, set_SD_rain, set_N_0rainD, set_n_moments, set_moments_file,
	#input
	profile_nlyr,profile_ngridx,profile_ngridy,nfreqs,freqs,
	profile_year,profile_month,profile_day,profile_time,
	profile_deltax,profile_deltay, profile_lat,profile_lon,profile_model_i,profile_model_j,
	profile_wind10u,profile_wind10v,profile_lfrac,
	profile_relhum_lev,profile_press_lev,profile_temp_lev,profile_hgt_lev,
	profile_iwv,profile_cwp,profile_iwp,profile_rwp,profile_swp,profile_gwp,
	profile_cwc_q,profile_iwc_q,profile_rwc_q,profile_swc_q,profile_gwc_q)
	
	return result_pamtraVersion, result_pamtraHash, result_Ze, result_attenuationHydro, result_attenuationAtmo, result_hgt, result_tb, result_angles


class pyPamtra:
	
	def __init__(self):
		
		self.set = dict()
		self.set["verbose"]=0
		self.set["dump_to_file"]=False
		self.set["tmp_path"]='/tmp/'
		self.set["data_path"]='/home/mech/models/pamtra/data/'

		self.set["obs_height"]=833000.
		self.set["units"]='T'
		self.set["outpol"]='VH'
		self.set["creator"]='Pamtrauser'

		self.set["active"]=True
		self.set["passive"]=True

		self.set["ground_type"]='S'
		self.set["salinity"]=33.0
		self.set["emissivity"]=0.6

		self.set["lgas_extinction"]=True
		self.set["gas_mod"]='R98'

		self.set["lhyd_extinction"]=True
		self.set["lphase_flag"]= True

		self.set["SD_snow"]='Exp' 
		self.set["N_0snowDsnow"]=7.628 
		self.set["EM_snow"]='icesf' 
		self.set["SP"]=0.2 
		self.set["isnow_n0"]=1
		self.set["liu_type"]=8

		self.set["SD_grau"]='Exp' 
		self.set["N_0grauDgrau"]=4.0 
		self.set["EM_grau"]='surus'

		self.set["EM_ice"]='mieic'

		self.set["SD_rain"]='Exp' 
		self.set["N_0rainD"]=8.0

		self.set["n_moments"]=1
		self.set["moments_file"]='snowCRYSTAL'
		
		self.setDefaultKeys = self.set.keys()
		
	
	def readProfile(self,inputFile):
		
		self.p = dict()
		
		f = open(inputFile,"r")
		g = csv.reader(f,delimiter=" ",skipinitialspace=True)
		self.p["year"], self.p["month"], self.p["day"], self.p["time"], self.p["ngridx"], self.p["ngridy"], self.p["nlyr"], self.p["deltax"], self.p["deltay"] = g.next()
		
		self.p["ngridx"] = int(self.p["ngridx"])
		self.p["ngridy"] = int(self.p["ngridy"])
		self.p["nlyr"] = int(self.p["nlyr"])
		
		shape2D = (self.p["ngridx"],self.p["ngridy"],)
		shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyr"],)
		shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyr"]+1,)
		
		self.p["model_i"] = np.zeros(shape2D)
		self.p["model_j"] = np.zeros(shape2D)
		self.p["lat"] = np.zeros(shape2D)
		self.p["lon"] = np.zeros(shape2D)
		self.p["lfrac"] = np.zeros(shape2D)
		self.p["wind10u"] = np.zeros(shape2D)
		self.p["wind10v"] = np.zeros(shape2D)
		
		self.p["iwv"] = np.zeros(shape2D)
		self.p["cwp"] = np.zeros(shape2D)
		self.p["iwp"] = np.zeros(shape2D)
		self.p["rwp"] = np.zeros(shape2D)
		self.p["swp"] = np.zeros(shape2D)
		self.p["gwp"] = np.zeros(shape2D)
		
		self.p["hgt_lev"] = np.zeros(shape3Dplus)
		self.p["temp_lev"] = np.zeros(shape3Dplus)
		self.p["press_lev"] = np.zeros(shape3Dplus)
		self.p["relhum_lev"] = np.zeros(shape3Dplus)
		
		self.p["cwc_q"] = np.zeros(shape3D)
		self.p["iwc_q"] = np.zeros(shape3D)
		self.p["rwc_q"] = np.zeros(shape3D)
		self.p["swc_q"] = np.zeros(shape3D)
		self.p["gwc_q"] = np.zeros(shape3D)
		
		
		for x in np.arange(self.p["ngridx"]):
			for y in np.arange(self.p["ngridy"]):
				self.p["model_i"][x,y], self.p["model_j"][x,y] = np.array(g.next(),dtype=int)
				self.p["lat"][x,y], self.p["lon"][x,y], self.p["lfrac"][x,y],self.p["wind10u"][x,y],self.p["wind10v"][x,y]  = np.array(g.next(),dtype=float)
				self.p["iwv"][x,y],self.p["cwp"][x,y],self.p["iwp"][x,y],self.p["rwp"][x,y],self.p["swp"][x,y],self.p["gwp"][x,y] = np.array(g.next(),dtype=float)
				self.p["hgt_lev"][x,y,0],self.p["press_lev"][x,y,0],self.p["temp_lev"][x,y,0],self.p["relhum_lev"][x,y,0] = np.array(g.next(),dtype=float)
				for z in np.arange(self.p["nlyr"]):
					self.p["hgt_lev"][x,y,z+1],self.p["press_lev"][x,y,z+1],self.p["temp_lev"][x,y,z+1],self.p["relhum_lev"][x,y,z+1],self.p["cwc_q"][x,y,z],self.p["iwc_q"][x,y,z],self.p["rwc_q"][x,y,z],self.p["swc_q"][x,y,z],self.p["gwc_q"][x,y,z] = np.array(g.next(),dtype=float)

		f.close()

	
	
	def createProfile(self,timestamp,lat,lon,lfrac,wind10u,wind10v,
			iwv,cwp,iwp,rwp,swp,gwp,
			hgt_lev,press_lev,temp_lev,relhum_lev,
			cwc_q,iwc_q,rwc_q,swc_q,gwc_q):
		
		self.p = dict()
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
		self.p["nlyr"] = np.shape(hgt_lev)[-1]
		shape2D = (self.p["ngridx"],self.p["ngridy"],)
		shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyr"],)
		shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyr"]+1,)
		
		if (type(timestamp) == int) or (type(timestamp) == float) :
			timestamp = datetime.datetime.utcfromtimestamp(timestamp)
		self.p["year"] = timestamp.strftime("%Y")
		self.p["month"] = timestamp.strftime("%m")
		self.p["day"] = timestamp.strftime("%d")
		self.p["time"] = timestamp.strftime("%H%M")
		self.p["deltax"] = 0.
		self.p["deltay"] = 0.
		self.p["lat"] = lat.reshape(shape2D)
		self.p["lon"] = lon.reshape(shape2D)
		self.p["lfrac"] = lfrac.reshape(shape2D)
		self.p["model_i"] = np.where(hgt_lev.reshape(shape3D))[0]
		self.p["model_j"] = np.where(hgt_lev.reshape(shape3D))[1]
		self.p["wind10u"] = wind10u.reshape(shape2D)
		self.p["wind10v"] = wind10v.reshape(shape2D)

		self.p["iwv"] = iwv.reshape(shape2D)
		self.p["cwp"] = cwp.reshape(shape2D)
		self.p["iwp"] = iwp.reshape(shape2D)
		self.p["rwp"] = rwp.reshape(shape2D)
		self.p["swp"] = swp.reshape(shape2D)
		self.p["gwp"] = gwp.reshape(shape2D)
		
		self.p["hgt_lev"] = hgt_lev.reshape(shape3Dplus)
		self.p["temp_lev"] = temp_lev.reshape(shape3Dplus)
		self.p["press_lev"] = press_lev.reshape(shape3Dplus)
		self.p["relhum_lev"] = relhum_lev.reshape(shape3Dplus)
		
		self.p["cwc_q"] = cwc_q.reshape(shape3D)
		self.p["iwc_q"] = iwc_q.reshape(shape3D)
		self.p["rwc_q"] = rwc_q.reshape(shape3D)
		self.p["swc_q"] = swc_q.reshape(shape3D)
		self.p["gwc_q"] = gwc_q.reshape(shape3D)

	def runParallelPamtra(self,freqs):
		import pp

		
		ppservers=("localhost",)
		job_server = pp.Server(ppservers=ppservers) 
		print "Starting pp with", job_server.get_ncpus(), "workers"
		
		

		if (type(freqs) == int) or (type(freqs) == float): freqs = [freqs]
		
		self.freqs = freqs
		self.nfreqs = len(freqs)
		
		for key in self.set:
			if key not in self.setDefaultKeys:
				print "Warning Could not parse ",key
				
		if self.set["n_moments"]==2:
			raise IOError("multi moments not implemented yet!")

		self.r = dict()

		job1 = job_server.submit(PamtraFortranWrapper, (
		#self.set
		self.set["verbose"], self.set["dump_to_file"], self.set["tmp_path"], self.set["data_path"], self.set["obs_height"], self.set["units"], self.set["outpol"], self.set["creator"], self.set["active"], self.set["passive"], self.set["ground_type"], self.set["salinity"], self.set["emissivity"], self.set["lgas_extinction"], self.set["gas_mod"], self.set["lhyd_extinction"], self.set["lphase_flag"], self.set["SD_snow"], self.set["N_0snowDsnow"], self.set["EM_snow"], self.set["SP"], self.set["isnow_n0"], self.set["liu_type"], self.set["SD_grau"], self.set["N_0grauDgrau"], self.set["EM_grau"], self.set["EM_ice"], self.set["SD_rain"], self.set["N_0rainD"], self.set["n_moments"], self.set["moments_file"],
		#input
		self.p["nlyr"],self.p["ngridx"],self.p["ngridy"],self.nfreqs,self.freqs,
		self.p["year"],self.p["month"],self.p["day"],self.p["time"],
		self.p["deltax"],self.p["deltay"], self.p["lat"],self.p["lon"],self.p["model_i"],self.p["model_j"],
		self.p["wind10u"],self.p["wind10v"],self.p["lfrac"],
		self.p["relhum_lev"],self.p["press_lev"],self.p["temp_lev"],self.p["hgt_lev"],
		self.p["iwv"],self.p["cwp"],self.p["iwp"],self.p["rwp"],self.p["swp"],self.p["gwp"],
		self.p["cwc_q"],self.p["iwc_q"],self.p["rwc_q"],self.p["swc_q"],self.p["gwc_q"],), tuple(), ("pyPamtraLib","import numpy as np",))
		
		test = job1()
		
		print test
		
		self.r["pamtraVersion"],self.r["pamtraHash"],\
		self.r["Ze"],self.r["attenuationHydro"],self.r["attenuationAtmo"],self.r["hgt"], self.r["tb"], self.r["angles"] = test
		
		self.r["Ze_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["attenuationHydro_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["attenuationAtmo_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["tb_dimensions"] = ["gridx","gridy","outlevels","angles","frequency","stokes"]
		
		self.r["settings"] = self.set
		
		self.r["pamtraVersion"] = self.r["pamtraVersion"].strip()
		self.r["pamtraHash"] = self.r["pamtraHash"].strip()
		
		#for key in self.__dict__.keys():
			#print key
			#print self.__dict__[key]

	def runPamtra(self,freqs):

		if (type(freqs) == int) or (type(freqs) == float): freqs = [freqs]
		
		self.freqs = freqs
		self.nfreqs = len(freqs)
		
		for key in self.set:
			if key not in self.setDefaultKeys:
				print "Warning Could not parse ",key
				
		if self.set["n_moments"]==2:
			raise IOError("multi moments not implemented yet!")

		self.r = dict()

		#output
		self.r["pamtraVersion"],self.r["pamtraHash"],\
		self.r["Ze"],self.r["attenuationHydro"],self.r["attenuationAtmo"],self.r["hgt"], self.r["tb"], self.r["angles"] = \
		PamtraFortranWrapper(
		#self.set
		self.set["verbose"], self.set["dump_to_file"], self.set["tmp_path"], self.set["data_path"], self.set["obs_height"], self.set["units"], self.set["outpol"], self.set["creator"], self.set["active"], self.set["passive"], self.set["ground_type"], self.set["salinity"], self.set["emissivity"], self.set["lgas_extinction"], self.set["gas_mod"], self.set["lhyd_extinction"], self.set["lphase_flag"], self.set["SD_snow"], self.set["N_0snowDsnow"], self.set["EM_snow"], self.set["SP"], self.set["isnow_n0"], self.set["liu_type"], self.set["SD_grau"], self.set["N_0grauDgrau"], self.set["EM_grau"], self.set["EM_ice"], self.set["SD_rain"], self.set["N_0rainD"], self.set["n_moments"], self.set["moments_file"],
		#input
		self.p["nlyr"],self.p["ngridx"],self.p["ngridy"],self.nfreqs,self.freqs,
		self.p["year"],self.p["month"],self.p["day"],self.p["time"],
		self.p["deltax"],self.p["deltay"], self.p["lat"],self.p["lon"],self.p["model_i"],self.p["model_j"],
		self.p["wind10u"],self.p["wind10v"],self.p["lfrac"],
		self.p["relhum_lev"],self.p["press_lev"],self.p["temp_lev"],self.p["hgt_lev"],
		self.p["iwv"],self.p["cwp"],self.p["iwp"],self.p["rwp"],self.p["swp"],self.p["gwp"],
		self.p["cwc_q"],self.p["iwc_q"],self.p["rwc_q"],self.p["swc_q"],self.p["gwc_q"])
		
		self.r["Ze_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["attenuationHydro_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["attenuationAtmo_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["tb_dimensions"] = ["gridx","gridy","outlevels","angles","frequency","stokes"]
		
		self.r["settings"] = self.set
		
		self.r["pamtraVersion"] = self.r["pamtraVersion"].strip()
		self.r["pamtraHash"] = self.r["pamtraHash"].strip()
		
		#for key in self.__dict__.keys():
			#print key
			#print self.__dict__[key]



	def writeResultsToNumpy(self,fname):
		try: 
			self.r
			self.r["pamtraVersion"]
		except:
			raise IOError ("run runPamtra first!")
			
		f = open(fname, "w")
		pickle.dump(self.r, f)
		f.close()

	def loadResultsFromNumpy(self,fname):
		try: 
			f = open(fname, "r")
			self.r = pickle.load(f)
			f.close()
		except:
			raise IOError ("Could not read data")
		
	#calendar.timegm(date.timetuple())

	 #datetime.datetime.utcfromtimestamp(unix)
		
	def writeResultsToNetCDF(self,fname,ncForm="NETCDF3_CLASSIC"):
		import netCDF4
		try: 
			self.r
			self.r["pamtraVersion"]
		except:
			raise IOError ("run runPamtra first!")

		cdfFile = netCDF4.Dataset(fname,"w",ncForm)
		
		#write meta data
		cdfFile.history = "Created with pyPamtra (Version: "+self.r["pamtraVersion"]+", Git Hash: "+self.r["pamtraHash"]+") by "+self.set["creator"]+" (University of Cologne, IGMK) at " + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
		cdfFile.data_time = self.p["year"]+"/"+self.p["month"]+"/"+self.p["day"]+"-"+self.p["time"][0:2]+":"+self.p["time"][2:4]

		#make dimesnions
		
		cdfFile.createDimension('nlon',self.p["ngridx"])
		cdfFile.createDimension('nlat',self.p["ngridy"])
		cdfFile.createDimension('nfreq',self.nfreqs)
		if (self.r["settings"]["passive"]):
			cdfFile.createDimension('nang',len(self.r["angles"]))
			cdfFile.createDimension('nout',np.shape(self.r["tb"])[2])
			cdfFile.createDimension('nstokes',np.shape(self.r["tb"])[-1])
		if (self.r["settings"]["active"]):
			cdfFile.createDimension('nlyr',self.p["nlyr"])
		
		#create variables
		nc_angle = cdfFile.createVariable('angle','f4',('nang',),fill_value= -9999)
		nc_angle.units = 'deg'
		
		nc_frequency = cdfFile.createVariable('frequency','f4',('nfreq',),fill_value= -9999)
		nc_frequency.units = 'GHz'

		dim2d = ("nlon","nlat",)
		nc_model_i = cdfFile.createVariable('model_i', 'i4',dim2d,fill_value= -9999)
		nc_model_i.units = "-"

		nc_model_j = cdfFile.createVariable('model_j', 'i4',dim2d,fill_value= -9999)
		nc_model_j.units = "-"

		nc_longitude = cdfFile.createVariable('longitude', 'f4',dim2d,fill_value= -9999)
		nc_longitude.units = "deg.dec"

		nc_latitude = cdfFile.createVariable('latitude', 'f4',dim2d,fill_value= -9999)
		nc_latitude.units = "deg.dec"

		nc_lfrac = cdfFile.createVariable('lfrac', 'f4',dim2d,fill_value= -9999)
		nc_lfrac.units = "-"

		nc_iwv = cdfFile.createVariable('iwv', 'f4',dim2d,fill_value= -9999)
		nc_iwv.units = "kg/m^2"

		nc_cwp= cdfFile.createVariable('cwp', 'f4',dim2d,fill_value= -9999)
		nc_cwp.units = "kg/m^2"

		nc_iwp = cdfFile.createVariable('iwp', 'f4',dim2d,fill_value= -9999)
		nc_iwp.units = "kg/m^2"

		nc_rwp = cdfFile.createVariable('rwp', 'f4',dim2d,fill_value= -9999)
		nc_rwp.units = "kg/m^2"

		nc_swp = cdfFile.createVariable('swp', 'f4',dim2d,fill_value= -9999)
		nc_swp.units = "kg/m^2"

		nc_gwp = cdfFile.createVariable('gwp', 'f4',dim2d,fill_value= -9999)
		nc_gwp.units = "kg/m^2"

		nc_hwp = cdfFile.createVariable('hwp', 'f4',dim2d,fill_value= -9999)
		nc_hwp.units = "kg/m^2"

		if (self.r["settings"]["active"]):

			dim3d = ("nlon","nlat","nlyr",)
			nc_height = cdfFile.createVariable('height', 'f4',dim3d,fill_value= -9999)
			nc_height.units = "m"

			dim4d = ("nlon","nlat","nlyr","nfreq")
			nc_Ze = cdfFile.createVariable('Ze', 'f4',dim4d,fill_value= -9999)
			nc_Ze.units = "dBz"

			nc_Attenuation_Hydrometeors = cdfFile.createVariable('Attenuation_Hydrometeors', 'f4',dim4d,fill_value= -9999)
			nc_Attenuation_Hydrometeors.units = "dB"

			nc_Attenuation_Atmosphere = cdfFile.createVariable('Attenuation_Atmosphere', 'f4',dim4d,fill_value= -9999)
			nc_Attenuation_Atmosphere.units = "dB"


		if (self.r["settings"]["passive"]):
			dim6d = ("nlon","nlat","nout","nang","nfreq","nstokes")
			nc_tb = cdfFile.createVariable('tb', 'f4',dim6d,fill_value= -9999)
			nc_tb.units = "K"

		#save data
		nc_angle[:] = self.r["angles"]
		nc_frequency[:] = self.freqs
		nc_model_i[:] = self.p["model_i"]
		nc_model_j[:] = self.p["model_j"]
		nc_longitude[:] = self.p["lon"]
		nc_latitude[:] = self.p["lat"]
		nc_lfrac[:] = self.p["lfrac"]
		nc_iwv[:] = self.p["iwv"]
		nc_cwp[:] = self.p["cwp"]
		nc_iwp[:] = self.p["iwp"]
		nc_rwp[:] = self.p["rwp"]
		nc_swp[:] = self.p["swp"]
		nc_gwp[:] = self.p["gwp"]
		if (self.r["settings"]["n_moments"] == 2): nc_hwp[:] = self.p["hwp"]
		
		if (self.r["settings"]["passive"]):
			nc_tb[:] = self.r["tb"]
		if (self.r["settings"]["active"]):          
			nc_height[:] = self.r["hgt"]
			nc_Ze[:] = self.r["Ze"]
			nc_Attenuation_Hydrometeors[:] = self.r["attenuationHydro"]
			nc_Attenuation_Atmosphere[:] = self.r["attenuationAtmo"]
			
		cdfFile.close()
		
		

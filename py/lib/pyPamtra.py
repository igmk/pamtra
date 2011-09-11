# -*- coding: utf-8 -*-

#http://www.parallelpython.com/component/option,com_smf/Itemid,29/topic,407.0

import numpy as np
import datetime
import csv
import pickle
import time,calendar, datetime
import warnings

import meteoSI

try:
	import pp
except:
	warnings.warn("parallel python not available", Warning)

missingNumber=-9999

def PamtraFortranWrapper(*args):
	#is needed because pp cannot work with fortran modules directly
	import pyPamtraLib
	code = ""
	for ii in range(len(args)):
		code = code + "args["+str(ii)+"],"
	return eval("pyPamtraLib.pypamtralib("+code[0:-1]+")")


class pyPamtra:
	
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
		
		self.dimensions = dict()
		
		self.dimensions["year"] = []
		self.dimensions["month"] = []
		self.dimensions["day"] = []
		self.dimensions["time"] = []
		
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
		
		self.dimensions["hgt_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
		self.dimensions["temp_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
		self.dimensions["p_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
		self.dimensions["relhum_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
		
		self.dimensions["cwc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["iwc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["rwc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["swc_q"] = ["ngridx","ngridy","max_nlyrs"]
		self.dimensions["gwc_q"] = ["ngridx","ngridy","max_nlyrs"]

		
		self.units = dict()
		
		self.units["timestamp"] = "seconds since 1970-01-01"
		self.units["year"] = "yyyy"
		self.units["month"] = "mm"
		self.units["day"] = "dd"
		self.units["time"] = "HHMM"
		
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
		
		self.units["hgt_lev"] = "m"
		self.units["temp_lev"] = "K"
		self.units["p_lev"] = "Pa"
		self.units["relhum_lev"] = "1"
		
		self.units["cwc_q"] = "kg/kg"
		self.units["iwc_q"] = "kg/kg"
		self.units["rwc_q"] = "kg/kg"
		self.units["swc_q"] = "kg/kg"
		self.units["gwc_q"] = "kg/kg"

		
		
		
	
	def readPamtraProfile(self,inputFile):
		
		self.p = dict()
		
		f = open(inputFile,"r")
		g = csv.reader(f,delimiter=" ",skipinitialspace=True)
		self.p["year"], self.p["month"], self.p["day"], self.p["time"], self.p["ngridx"], self.p["ngridy"], self.p["nlyrs"], self.p["deltax"], self.p["deltay"] = g.next()
		
		self.p["ngridx"] = int(self.p["ngridx"])
		self.p["ngridy"] = int(self.p["ngridy"])
		self.p["nlyrs"] = int(self.p["nlyrs"])
		
		shape2D = (self.p["ngridx"],self.p["ngridy"],)
		shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],)
		shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"]+1,)
		
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
				for z in np.arange(self.p["nlyrs"]):
					self.p["hgt_lev"][x,y,z+1],self.p["press_lev"][x,y,z+1],self.p["temp_lev"][x,y,z+1],self.p["relhum_lev"][x,y,z+1],self.p["cwc_q"][x,y,z],self.p["iwc_q"][x,y,z],self.p["rwc_q"][x,y,z],self.p["swc_q"][x,y,z],self.p["gwc_q"][x,y,z] = np.array(g.next(),dtype=float)

		#in PyPamtra I want relhum not in %
		self.p["relhum_lev"] = self.p["relhum_lev"]/100.

		#finally make an array from nlyrs
		self.p["nlyrs"] = np.ones(shape2D,dtype=int)*self.p["nlyrs"]

		f.close()


	def createProfile(self,timestamp,lat,lon,lfrac,wind10u,wind10v,
			hgt_lev,press_lev,temp_lev,relhum_lev,
			cwc_q,iwc_q,rwc_q,swc_q,gwc_q):
		
		q_lev = meteoSI.rh2q(relhum_lev,temp_lev,press_lev)
		q = (q_lev[...,0:-1] + q_lev[...,1:])/2.
		print "can I average p linearly?"
		rho_moist_lev = meteoSI.moist_rho_q(press_lev,temp_lev,q_lev)
		rho_moist = (rho_moist_lev[...,0:-1] + rho_moist_lev[...,1:])/2.

		dz = np.diff(hgt_lev)

		iwv = meteoSI.integrate_xq_xwp(q,rho_moist,hgt_lev)
		cwp = meteoSI.integrate_xq_xwp(cwc_q,rho_moist,hgt_lev)
		iwp = meteoSI.integrate_xq_xwp(iwc_q,rho_moist,hgt_lev)
		rwp = meteoSI.integrate_xq_xwp(rwc_q,rho_moist,hgt_lev)
		swp = meteoSI.integrate_xq_xwp(swc_q,rho_moist,hgt_lev)
		gwp = meteoSI.integrate_xq_xwp(gwc_q,rho_moist,hgt_lev)
		
		
		return self.createFullProfile(timestamp,lat,lon,lfrac,wind10u,wind10v,
		iwv,cwp,iwp,rwp,swp,gwp,
		hgt_lev,press_lev,temp_lev,relhum_lev,
		cwc_q,iwc_q,rwc_q,swc_q,gwc_q,missingNumber=missingNumber)



	
	def createFullProfile(self,timestamp,lat,lon,lfrac,wind10u,wind10v,
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
		self.p["nlyr"] = np.shape(hgt_lev)[-1] -1
		shape2D = (self.p["ngridx"],self.p["ngridy"],)
		shape3D = (self.p["ngridx"],self.p["ngridy"],np.max(self.p["nlyrs"]),)
		shape3Dplus = (self.p["ngridx"],self.p["ngridy"],p.max(self.p["nlyrs"])+1,)
		
		
		if (type(timestamp) == str):
			self.p["year"] = timestamp[0:4]
			self.p["month"] = timestamp[4:6]
			self.p["day"] = timestamp[6:8]
			self.p["time"] = timestamp[8:12]
		else:
			if (type(timestamp) == int) or (type(timestamp) == float) :
				timestamp = datetime.datetime.utcfromtimestamp(timestamp)
			elif (type(timestamp) == datetime):
				pass
			else:
				raise TypeError("tiemstamp has to be int or float or string")
			self.p["year"] = timestamp.strftime("%Y")
			self.p["month"] = timestamp.strftime("%m")
			self.p["day"] = timestamp.strftime("%d")
			self.p["time"] = timestamp.strftime("%H%M")
				
		self.p["deltax"] = 0.
		self.p["deltay"] = 0.
		self.p["lat"] = lat.reshape(shape2D)
		self.p["lon"] = lon.reshape(shape2D)
		self.p["lfrac"] = lfrac.reshape(shape2D)
		self.p["model_i"] = np.array(np.where(lat.reshape(shape2D))[0]).reshape(shape2D) +1
		self.p["model_j"] = np.array(np.where(lon.reshape(shape2D))[1]).reshape(shape2D) +1
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

	def runParallelPamtra(self,freqs,pp_servers=(),pp_local_workers=1,pp_deltaF=0,pp_deltaX=0,pp_deltaY = 0):
		
		if pp_local_workers == "auto":
			job_server = pp.Server(ppservers=pp_servers) 
		else:
			job_server = pp.Server(pp_local_workers,ppservers=pp_servers) 
		if self.set["pyVerbose"] > = 0: 
			print "Starting pp with: "
			pp_nodes = job_server.get_active_nodes()
			for key in pp_nodes.keys():
				print key+": "+str(pp_nodes[key])+" nodes"
		
		if (type(freqs) == int) or (type(freqs) == float): freqs = [freqs]
		
		self.freqs = freqs
		self.nfreqs = len(freqs)
		
		self.nstokes = 2
		self.noutlevels = 2
		self.nangles = 32

		
		for key in self.set:
			if key not in self.setDefaultKeys:
				warnings.warn("Warning can not parse setting: ",key, Warning)
				
		if self.set["n_moments"]==2:
			raise IOError("multi moments not implemented yet!")

		self.r = dict()

		
		if pp_deltaF==0: pp_deltaF = self.nfreqs
		if pp_deltaX==0: pp_deltaX = self.p["ngridx"]
		if pp_deltaY==0: pp_deltaY = self.p["ngridy"]
		
		pp_ii = -1
		pp_jobs = dict()
		
		self.r["Ze"] = np.ones((self.p["ngridx"],self.p["ngridy"],np.max(self.p["nlyrs"]),self.nfreqs,))*missingNumber
		self.r["attenuationHydro"] = np.ones((self.p["ngridx"],self.p["ngridy"],np.max(self.p["nlyrs"]),self.nfreqs,))*missingNumber
		self.r["attenuationAtmo"] = np.ones((self.p["ngridx"],self.p["ngridy"],np.max(self.p["nlyrs"]),self.nfreqs,))*missingNumber
		self.r["hgt"] = np.ones((self.p["ngridx"],self.p["ngridy"],np.max(self.p["nlyrs"]),))*missingNumber
		self.r["tb"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.noutlevels,self.nangles,self.nfreqs,self.nstokes))*missingNumber
		
		self.r["Ze_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["attenuationHydro_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["attenuationAtmo_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["hgt_dimensions"] = ["gridx","gridy","lyr"]
		self.r["tb_dimensions"] = ["gridx","gridy","outlevels","angles","frequency","stokes"]

		
		for pp_startF in np.arange(0,self.nfreqs,pp_deltaF):
			pp_endF = pp_startF + pp_deltaF
			if pp_endF > self.nfreqs: pp_endF = self.nfreqs
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
					
					pp_jobs[pp_ii] = job_server.submit(PamtraFortranWrapper, (
					#self.set
					self.set["verbose"], self.set["dump_to_file"], self.set["tmp_path"], self.set["data_path"], self.set["obs_height"], self.set["units"], self.set["outpol"], self.set["creator"], self.set["active"], self.set["passive"], self.set["ground_type"], self.set["salinity"], self.set["emissivity"], self.set["lgas_extinction"], self.set["gas_mod"], self.set["lhyd_extinction"], self.set["lphase_flag"], self.set["SD_snow"], self.set["N_0snowDsnow"], self.set["EM_snow"], self.set["SP"], self.set["isnow_n0"], self.set["liu_type"], self.set["SD_grau"], self.set["N_0grauDgrau"], self.set["EM_grau"], self.set["EM_ice"], self.set["SD_rain"], self.set["N_0rainD"], self.set["n_moments"], self.set["moments_file"],
					#input
					pp_ngridx,
					pp_ngridy,
					np.max(self.p["nlyrs"]),
					self.p["nlyrs"][pp_startX:pp_endX,pp_startY:pp_endY],
					pp_nfreqs,
					self.freqs[pp_startF:pp_endF],
					self.p["year"],self.p["month"],self.p["day"],self.p["time"],
					self.p["deltax"],self.p["deltay"],
					self.p["lat"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["lon"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["model_i"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["model_j"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["wind10u"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["wind10v"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["lfrac"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["relhum_lev"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["press_lev"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["temp_lev"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["hgt_lev"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["iwv"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["cwp"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["iwp"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["rwp"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["swp"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["gwp"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["cwc_q"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["iwc_q"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["rwc_q"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["swc_q"][pp_startX:pp_endX,pp_startY:pp_endY],
					self.p["gwc_q"][pp_startX:pp_endX,pp_startY:pp_endY],
					),tuple(), ("pyPamtraLib","numpy",))
										
		job_server.wait()
		pp_ii = -1

		for pp_startF in np.arange(0,self.nfreqs,pp_deltaF):
			pp_endF = pp_startF + pp_deltaF
			if pp_endF > self.nfreqs: pp_endF = self.nfreqs
			for pp_startX in np.arange(0,self.p["ngridx"],pp_deltaX):
				pp_endX = pp_startX + pp_deltaX
				if pp_endX > self.p["ngridx"]: pp_endX = self.p["ngridx"]
				for pp_startY in np.arange(0,self.p["ngridy"],pp_deltaY):
					pp_endY = pp_startY + pp_deltaY
					if pp_endY > self.p["ngridy"]: pp_endY = self.p["ngridy"]
					pp_ii +=1

					(self.r["pamtraVersion"],self.r["pamtraHash"],
					self.r["Ze"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF],
					self.r["attenuationHydro"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF],
					self.r["attenuationAtmo"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF],
					self.r["hgt"][pp_startX:pp_endX,pp_startY:pp_endY,:], 
					self.r["tb"][pp_startX:pp_endX,pp_startY:pp_endY,:,:,pp_startF:pp_endF,:], 
					self.r["angles"] ) = pp_jobs[pp_ii]()
		

		
		self.r["settings"] = self.set
		
		self.r["pamtraVersion"] = self.r["pamtraVersion"].strip()
		self.r["pamtraHash"] = self.r["pamtraHash"].strip()
		
		#for key in self.__dict__.keys():
			#print key
			#print self.__dict__[key]
		if self.set["pyVerbose"] > = 0: job_server.print_stats()
		job_server.destroy()
		
	def runPamtra(self,freqs):
		
		if (type(freqs) == int) or (type(freqs) == float): freqs = [freqs]
		
		self.freqs = freqs
		self.nfreqs = len(freqs)
		
		for key in self.set:
			if key not in self.setDefaultKeys:
				print "Warning Could not parse ",key
				
		if self.set["n_moments"]==2:
			raise IOError("multi moments not implemented yet!")
		
		if np.max(self.p["relhum_lev"])>1.5:
			raise IOError("relative humidity is _not_ in %!")
		
		tttt = time.time()
		self.r = dict()

		#output
		self.r["pamtraVersion"],self.r["pamtraHash"],\
		self.r["Ze"],self.r["attenuationHydro"],self.r["attenuationAtmo"],self.r["hgt"], self.r["tb"], self.r["angles"] = \
		PamtraFortranWrapper(
		#self.set
		self.set["verbose"], self.set["dump_to_file"], self.set["tmp_path"], self.set["data_path"], self.set["obs_height"], self.set["units"], self.set["outpol"], self.set["creator"], self.set["active"], self.set["passive"], self.set["ground_type"], self.set["salinity"], self.set["emissivity"], self.set["lgas_extinction"], self.set["gas_mod"], self.set["lhyd_extinction"], self.set["lphase_flag"], self.set["SD_snow"], self.set["N_0snowDsnow"], self.set["EM_snow"], self.set["SP"], self.set["isnow_n0"], self.set["liu_type"], self.set["SD_grau"], self.set["N_0grauDgrau"], self.set["EM_grau"], self.set["EM_ice"], self.set["SD_rain"], self.set["N_0rainD"], self.set["n_moments"], self.set["moments_file"],
		#input
		self.p["ngridx"],self.p["ngridy"],np.max(self.p["nlyrs"]),self.p["nlyrs"],self.nfreqs,self.freqs,
		self.p["year"],self.p["month"],self.p["day"],self.p["time"],
		self.p["deltax"],self.p["deltay"], self.p["lat"],self.p["lon"],self.p["model_i"],self.p["model_j"],
		self.p["wind10u"],self.p["wind10v"],self.p["lfrac"],
		self.p["relhum_lev"],self.p["press_lev"],self.p["temp_lev"],self.p["hgt_lev"],
		self.p["iwv"],self.p["cwp"],self.p["iwp"],self.p["rwp"],self.p["swp"],self.p["gwp"],
		self.p["cwc_q"],self.p["iwc_q"],self.p["rwc_q"],self.p["swc_q"],self.p["gwc_q"])
		
		self.dimensions["Ze"] = ["gridx","gridy","lyr","frequency"]
		self.dimensions["attenuationHydros"] = ["gridx","gridy","lyr","frequency"]
		self.dimensions["attenuationAtmo"] = ["gridx","gridy","lyr","frequency"]
		self.dimensions["tb"] = ["gridx","gridy","outlevels","angles","frequency","stokes"]
		
		self.r["settings"] = self.set
		
		self.r["pamtraVersion"] = self.r["pamtraVersion"].strip()
		self.r["pamtraHash"] = self.r["pamtraHash"].strip()
		
		#for key in self.__dict__.keys():
			#print key
			#print self.__dict__[key]
		
		if self.set["pyVerbose"] > = 0: print "pyPamtra runtime:", time.time() - tttt



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
			cdfFile.createDimension('nlyr',np.max(self.p["nlyrs"]))
		
		#create variables
		nc_angle = cdfFile.createVariable('angle','f4',('nang',),fill_value= missingNumber)
		nc_angle.units = 'deg'
		
		nc_frequency = cdfFile.createVariable('frequency','f4',('nfreq',),fill_value= missingNumber)
		nc_frequency.units = 'GHz'

		dim2d = ("nlon","nlat",)
		nc_model_i = cdfFile.createVariable('model_i', 'i4',dim2d,fill_value= missingNumber)
		nc_model_i.units = "-"

		nc_model_j = cdfFile.createVariable('model_j', 'i4',dim2d,fill_value= missingNumber)
		nc_model_j.units = "-"

		nc_nlyrs = cdfFile.createVariable('nlyrs', 'i4',dim2d,fill_value= missingNumber)
		nc_nlyrs.units = "-"

		nc_longitude = cdfFile.createVariable('longitude', 'f4',dim2d,fill_value= missingNumber)
		nc_longitude.units = "deg.dec"

		nc_latitude = cdfFile.createVariable('latitude', 'f4',dim2d,fill_value= missingNumber)
		nc_latitude.units = "deg.dec"

		nc_lfrac = cdfFile.createVariable('lfrac', 'f4',dim2d,fill_value= missingNumber)
		nc_lfrac.units = "-"

		nc_iwv = cdfFile.createVariable('iwv', 'f4',dim2d,fill_value= missingNumber)
		nc_iwv.units = "kg/m^2"

		nc_cwp= cdfFile.createVariable('cwp', 'f4',dim2d,fill_value= missingNumber)
		nc_cwp.units = "kg/m^2"

		nc_iwp = cdfFile.createVariable('iwp', 'f4',dim2d,fill_value= missingNumber)
		nc_iwp.units = "kg/m^2"

		nc_rwp = cdfFile.createVariable('rwp', 'f4',dim2d,fill_value= missingNumber)
		nc_rwp.units = "kg/m^2"

		nc_swp = cdfFile.createVariable('swp', 'f4',dim2d,fill_value= missingNumber)
		nc_swp.units = "kg/m^2"

		nc_gwp = cdfFile.createVariable('gwp', 'f4',dim2d,fill_value= missingNumber)
		nc_gwp.units = "kg/m^2"

		#nc_hwp = cdfFile.createVariable('hwp', 'f4',dim2d,fill_value= missingNumber)
		#nc_hwp.units = "kg/m^2"

		if (self.r["settings"]["active"]):

			dim3d = ("nlon","nlat","nlyr",)
			nc_height = cdfFile.createVariable('height', 'f4',dim3d,fill_value= missingNumber)
			nc_height.units = "m"

			dim4d = ("nlon","nlat","nlyr","nfreq")
			nc_Ze = cdfFile.createVariable('Ze', 'f4',dim4d,fill_value= missingNumber)
			nc_Ze.units = "dBz"

			nc_Attenuation_Hydrometeors = cdfFile.createVariable('Attenuation_Hydrometeors', 'f4',dim4d,fill_value= missingNumber)
			nc_Attenuation_Hydrometeors.units = "dB"

			nc_Attenuation_Atmosphere = cdfFile.createVariable('Attenuation_Atmosphere', 'f4',dim4d,fill_value= missingNumber)
			nc_Attenuation_Atmosphere.units = "dB"


		if (self.r["settings"]["passive"]):
			dim6d = ("nlon","nlat","nout","nang","nfreq","nstokes")
			nc_tb = cdfFile.createVariable('tb', 'f4',dim6d,fill_value= missingNumber)
			nc_tb.units = "K"

		#save data
		nc_angle[:] = self.r["angles"]
		nc_frequency[:] = self.freqs
		nc_nlyrs[:] = self.p["nlyrs"]
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
		#if (self.r["settings"]["n_moments"] == 2): nc_hwp[:] = self.p["hwp"]
		
		if (self.r["settings"]["passive"]):
			nc_tb[:] = self.r["tb"]
		if (self.r["settings"]["active"]):
			nc_height[:] = self.r["hgt"]
			nc_Ze[:] = self.r["Ze"]
			nc_Attenuation_Hydrometeors[:] = self.r["attenuationHydro"]
			nc_Attenuation_Atmosphere[:] = self.r["attenuationAtmo"]
			
		cdfFile.close()



# -*- coding: utf-8 -*-

#import numpy as np

import pyPamtraLib


class pyPamtra:
		#def __init__(self,z,p,t,rh,qv,cwc,iwc,rwc,swc,gwc,year,month,day,time,lat,lon):
		##print yyyy, mm, dd, hhmm, nx, ny, nz, dx, dy
		##print lat, lon, lfrac, wind10
		#self.levels = len(z[0,:])
		#self.nprof =  len(z[:,0])
		#shape = (self.nprof,self.levels)
		#shape1 = (self.nprof,self.levels-1)
		#self.year = year
		#self.mon = month
		#self.day = day
		#self.time = time
		#self.nx = self.nprof
		#self.ny = 1
		#self.nz = self.levels - 1
		#self.dx = 0.
		#self.dy = 0.
		#self.lat = lat
		#self.lon = lon
		#self.lfrac = -1.0
		#self.windu = 0.0
		#self.windv = 0.0
		#self.z = z.reshape(shape).astype(float)
		#self.p = p.reshape(shape).astype(float)
		#self.t = t.reshape(shape).astype(float)
		#self.rh = rh.reshape(shape).astype(float)
		#self.qv = qv.reshape(shape).astype(float)
		#self.cwc = cwc.reshape(shape1).astype(float)
		#self.iwc = iwc.reshape(shape1).astype(float)
		#self.rwc = rwc.reshape(shape1).astype(float)
		#self.swc = swc.reshape(shape1).astype(float)
		#self.gwc = gwc.reshape(shape1).astype(float)
		## Calculate additional values
		#self.cwc_q = np.zeros(shape1)
		#self.iwc_q = np.zeros(shape1)
		#self.rwc_q = np.zeros(shape1)
		#self.swc_q = np.zeros(shape1)
		#self.gwc_q = np.zeros(shape1)
		#self.delz = np.zeros(shape1)
		#self.wvc = np.zeros(shape1)
		
		#for i in range(self.levels-1):
		## end for
	## end def __init__
	#for j in range(self.nprof):
		#self.delz[j,i] = self.z[j,i+1]-self.z[j,i]
		#zbar = (self.z[j,i+1]+self.z[j,i])/2.
		#tbar = (self.t[j,i+1]+self.t[j,i])/2.
		#rhbar = (self.rh[j,i+1]+self.rh[j,i])/2.
		#xp = -1.*np.log(self.p[j,i+1]/self.p[j,i])/self.delz[j,i]
		#pbar = -1.*self.p[j,i]/xp*(np.exp(-xp*self.delz[j,i])-1.)/self.delz[j,i]
		#rho_m = moist_rho(pbar,tbar,rh2q(rhbar,tbar,pbar))
		##import pdb; pdb.set_trace()
		#self.cwc_q[j,i] = self.cwc[j,i]/rho_m
		#self.iwc_q[j,i] = self.iwc[j,i]/rho_m
		#self.rwc_q[j,i] = self.rwc[j,i]/rho_m
		#self.swc_q[j,i] = self.swc[j,i]/rho_m
		#self.gwc_q[j,i] = self.gwc[j,i]/rho_m
		#self.wvc[j,i] = (self.qv[j,i+1]+self.qv[j,i])/2.*rho_m

	
	def runPamtra(self,inputFile,frequencies,userSettings):

		if (type(freqs) == int) or ((type(freqs) == float): freqs = [freqs]

		settings = dict()
		settings["verbose"]=0
		settings["write_nc"]=True
		settings["dump_to_file"]=False
		settings["tmp_path"]='/tmp/'
		settings["write_nc"]=True
		settings["input_path"]='../test/referenceProfile'
		settings["output_path"]='../test/tmp'
		settings["data_path"]='/home/mech/models/pamtra/data/'

		settings["obs_height"]=833000.
		settings["units"]='T'
		settings["outpol"]='VH'
		settings["freq_str"]=''
		settings["file_desc"]=''
		settings["creator"]='Pamtrauser'

		settings["active"]=True
		settings["passive"]=True

		settings["ground_type"]='S'
		settings["salinity"]=33.0
		settings["emissivity"]=0.6

		settings["lgas_extinction"]=True
		settings["gas_mod"]='R98'

		settings["lhyd_extinction"]=True
		settings["lphase_flag"]= True

		settings["SD_snow"]='Exp' 
		settings["N_0snowDsnow"]=7.628 
		settings["EM_snow"]='icesf' 
		settings["SP"]=0.2 
		settings["isnow_n0"]=1
		settings["liu_type"]=8

		settings["SD_grau"]='Exp' 
		settings["N_0grauDgrau"]=4.0 
		settings["EM_grau"]='surus'

		settings["EM_ice"]='mieic'

		settings["SD_rain"]='Exp' 
		settings["N_0rainD"]=8.0

		settings["n_moments"]=1
		settings["moments_file"]='snowCRYSTAL'

		for key in userSettings:
			try: settings[key] = userSettings[key]
			except: "Could not parse ",key
		
		if settings["write_nc"] == False:
			raise ValueError("write_nc must be set true for historical reasons")
		
		
		nlyr = 50
		ngridx = 4
		ngridy = 1
		nfreqs = len(freqs)
		
		#output
		self.pamtraVersion,self.pamtraHash,\
		self.Ze,self.attenuationHydro,self.attenuationAtmo,self.hgt, self.tb, self.angles = \
		pyPamtraLib.pypamtralib(inputFile,
		#settings
		settings["verbose"], settings["write_nc"], settings["dump_to_file"], settings["input_path"], settings["output_path"], settings["tmp_path"], settings["data_path"], settings["obs_height"], settings["units"], settings["outpol"], settings["freq_str"], settings["file_desc"], settings["creator"], settings["active"], settings["passive"], settings["ground_type"], settings["salinity"], settings["emissivity"], settings["lgas_extinction"], settings["gas_mod"], settings["lhyd_extinction"], settings["lphase_flag"], settings["SD_snow"], settings["N_0snowDsnow"], settings["EM_snow"], settings["SP"], settings["isnow_n0"], settings["liu_type"], settings["SD_grau"], settings["N_0grauDgrau"], settings["EM_grau"], settings["EM_ice"], settings["SD_rain"], settings["N_0rainD"], settings["n_moments"], settings["moments_file"],
		#input
		nlyr,ngridx,ngridy,nfreqs,freqs)
		
		self.Ze_dimensions = ["gridx","gridy","lyr","frequency"]
		self.attenuationHydro_dimensions = ["gridx","gridy","lyr","frequency"]
		self.attenuationAtmo_dimensions = ["gridx","gridy","lyr","frequency"]
		self.tb_dimensions = ["gridx","gridy","outlevels","angles","frequency","stokes"]
		
		for key in self.__dict__.keys():
			print key
			print self.__dict__[key]
		
		
		

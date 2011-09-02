# -*- coding: utf-8 -*-

import numpy as np
import datetime
import pyPamtraLib
import csv
import pickle

class pyPamtra:
	
	def __init__(self):
		
		self.settings = dict()
		self.settings["verbose"]=0
		self.settings["dump_to_file"]=False
		self.settings["tmp_path"]='/tmp/'
		self.settings["data_path"]='/home/mech/models/pamtra/data/'

		self.settings["obs_height"]=833000.
		self.settings["units"]='T'
		self.settings["outpol"]='VH'
		self.settings["creator"]='Pamtrauser'

		self.settings["active"]=True
		self.settings["passive"]=True

		self.settings["ground_type"]='S'
		self.settings["salinity"]=33.0
		self.settings["emissivity"]=0.6

		self.settings["lgas_extinction"]=True
		self.settings["gas_mod"]='R98'

		self.settings["lhyd_extinction"]=True
		self.settings["lphase_flag"]= True

		self.settings["SD_snow"]='Exp' 
		self.settings["N_0snowDsnow"]=7.628 
		self.settings["EM_snow"]='icesf' 
		self.settings["SP"]=0.2 
		self.settings["isnow_n0"]=1
		self.settings["liu_type"]=8

		self.settings["SD_grau"]='Exp' 
		self.settings["N_0grauDgrau"]=4.0 
		self.settings["EM_grau"]='surus'

		self.settings["EM_ice"]='mieic'

		self.settings["SD_rain"]='Exp' 
		self.settings["N_0rainD"]=8.0

		self.settings["n_moments"]=1
		self.settings["moments_file"]='snowCRYSTAL'
		
		self.settingsDefaultKeys = self.settings.keys()
		
	
	def readProfile(self,inputFile):
		f = open(inputFile,"r")
		g = csv.reader(f,delimiter=" ",skipinitialspace=True)
		self.year, self.month, self.day, self.time, self.ngridx, self.ngridy, self.nlyr, self.deltax, self.deltay = g.next()
		
		self.ngridx = int(self.ngridx)
		self.ngridy = int(self.ngridy)
		self.nlyr = int(self.nlyr)
		
		shape2D = (self.ngridx,self.ngridy,)
		shape3D = (self.ngridx,self.ngridy,self.nlyr,)
		shape3Dplus = (self.ngridx,self.ngridy,self.nlyr+1,)
		
		self.model_i = np.zeros(shape2D)
		self.model_j = np.zeros(shape2D)
		self.lat = np.zeros(shape2D)
		self.lon = np.zeros(shape2D)
		self.lfrac = np.zeros(shape2D)
		self.wind10u = np.zeros(shape2D)
		self.wind10v = np.zeros(shape2D)
		
		self.iwv = np.zeros(shape2D)
		self.cwp = np.zeros(shape2D)
		self.iwp = np.zeros(shape2D)
		self.rwp = np.zeros(shape2D)
		self.swp = np.zeros(shape2D)
		self.gwp = np.zeros(shape2D)
		
		self.hgt_lev = np.zeros(shape3Dplus)
		self.temp_lev = np.zeros(shape3Dplus)
		self.press_lev = np.zeros(shape3Dplus)
		self.relhum_lev = np.zeros(shape3Dplus)
		
		self.cwc_q = np.zeros(shape3D)
		self.iwc_q = np.zeros(shape3D)
		self.rwc_q = np.zeros(shape3D)
		self.swc_q = np.zeros(shape3D)
		self.gwc_q = np.zeros(shape3D)
		
		
		for x in np.arange(self.ngridx):
			for y in np.arange(self.ngridy):
				print x,y
				self.model_i[x,y], self.model_j[x,y] = np.array(g.next(),dtype=int)
				self.lat[x,y], self.lon[x,y], self.lfrac[x,y],self.wind10u[x,y],self.wind10v[x,y]  = np.array(g.next(),dtype=float)
				self.iwv[x,y],self.cwp[x,y],self.iwp[x,y],self.rwp[x,y],self.swp[x,y],self.gwp[x,y] = np.array(g.next(),dtype=float)
				self.hgt_lev[x,y,0],self.press_lev[x,y,0],self.temp_lev[x,y,0],self.relhum_lev[x,y,0] = np.array(g.next(),dtype=float)
				for z in np.arange(self.nlyr):
					self.hgt_lev[x,y,z+1],self.press_lev[x,y,z+1],self.temp_lev[x,y,z+1],self.relhum_lev[x,y,z+1],self.cwc_q[x,y,z],self.iwc_q[x,y,z],self.rwc_q[x,y,z],self.swc_q[x,y,z],self.gwc_q[x,y,z] = np.array(g.next(),dtype=float)

		f.close()
		for key in self.__dict__.keys():
			print key,type(self.__dict__[key]),self.__dict__[key]

	
	
	def createProfile(self,timestamp,lat,lon,lfrac,wind10u,wind10v,
			iwv,cwp,iwp,rwp,swp,gwp,
			hgt_lev,press_lev,temp_lev,relhum_lev,
			cwc_q,iwc_q,rwc_q,swc_q,gwc_q):
		
		noDims = len(np.shape(hgt_lev))
		
		if noDims == 1:
			self.ngridx = 1
			self.ngridy = 1
		elif noDims == 2:
			self.ngridx = np.shape(hgt_lev)[0]
			self.ngridy = 1
		elif noDims == 3:
			self.ngridx = np.shape(hgt_lev)[0]
			self.ngridy = np.shape(hgt_lev)[1]
		else:
			print "Too many dimensions!"
			raise IOError
		self.nlyr = np.shape(hgt_lev)[-1]
		shape2D = (self.ngridx,self.ngridy,)
		shape3D = (self.ngridx,self.ngridy,self.nlyr,)
		shape3Dplus = (self.ngridx,self.ngridy,self.nlyr+1,)
		
		if (type(timestamp) == int) or (type(timestamp) == float) :
			timestamp = datetime.datetime.utcfromtimestamp(timestamp)
		self.year = timestamp.strftime("%Y")
		self.month = timestamp.strftime("%m")
		self.day = timestamp.strftime("%d")
		self.time = timestamp.strftime("%H%M")
		self.deltax = 0.
		self.deltay = 0.
		self.lat = lat.reshape(shape2D)
		self.lon = lon.reshape(shape2D)
		self.lfrac = lfrac.reshape(shape2D)
		self.model_i = np.where(hgt_lev.reshape(shape3D))[0]
		self.model_j = np.where(hgt_lev.reshape(shape3D))[1]
		self.wind10u = wind10u.reshape(shape2D)
		self.wind10v = wind10v.reshape(shape2D)

		self.iwv = iwv.reshape(shape2D)
		self.cwp = cwp.reshape(shape2D)
		self.iwp = iwp.reshape(shape2D)
		self.rwp = rwp.reshape(shape2D)
		self.swp = swp.reshape(shape2D)
		self.gwp = gwp.reshape(shape2D)
		
		self.hgt_lev = hgt_lev.reshape(shape3Dplus)
		self.temp_lev = temp_lev.reshape(shape3Dplus)
		self.press_lev = press_lev.reshape(shape3Dplus)
		self.relhum_lev = relhum_lev.reshape(shape3Dplus)
		
		self.cwc_q = cwc_q.reshape(shape3D)
		self.iwc_q = iwc_q.reshape(shape3D)
		self.rwc_q = rwc_q.reshape(shape3D)
		self.swc_q = swc_q.reshape(shape3D)
		self.gwc_q = gwc_q.reshape(shape3D)

	def runPamtra(self,freqs):

		if (type(freqs) == int) or (type(freqs) == float): freqs = [freqs]
		
		self.freqs = freqs
		self.nfreqs = len(freqs)
		
		for key in self.settings:
			if key not in self.settingsDefaultKeys:
				print "Warning Could not parse ",key

		#self.year = "2010"
		#self.month = "05"
		#self.day = "02"
		#self.time ="1200"

		#test = np.array([[ 1], [ 2], [ 3], [ 4]])

		#self.deltax=self.deltay= 2
		#self.lat=self.lon=self.model_i=self.model_j=test
		#self.wind10u=self.wind10v=self.lfrac=test*2
		#self.relhum_lev=self.press_lev=self.temp_lev=self.hgt_lev = np.random.random((self.ngridx,self.ngridy,self.nlyr+1))
		#self.iwv=self.cwp=self.iwp=self.rwp=self.swp=self.gwp = np.random.random((self.ngridx,self.ngridy))
		#self.cwc_q = self.iwc_q = self.rwc_q = self.swc_q = self.gwc_q = np.random.random((self.ngridx,self.ngridy,self.nlyr))
	
		self.r = dict()

		
		#output
		self.r["pamtraVersion"],self.r["pamtraHash"],\
		self.r["Ze"],self.r["attenuationHydro"],self.r["attenuationAtmo"],self.r["hgt"], self.r["tb"], self.r["angles"] = \
		pyPamtraLib.pypamtralib(
		#self.settings
		self.settings["verbose"], self.settings["dump_to_file"], self.settings["tmp_path"], self.settings["data_path"], self.settings["obs_height"], self.settings["units"], self.settings["outpol"], self.settings["creator"], self.settings["active"], self.settings["passive"], self.settings["ground_type"], self.settings["salinity"], self.settings["emissivity"], self.settings["lgas_extinction"], self.settings["gas_mod"], self.settings["lhyd_extinction"], self.settings["lphase_flag"], self.settings["SD_snow"], self.settings["N_0snowDsnow"], self.settings["EM_snow"], self.settings["SP"], self.settings["isnow_n0"], self.settings["liu_type"], self.settings["SD_grau"], self.settings["N_0grauDgrau"], self.settings["EM_grau"], self.settings["EM_ice"], self.settings["SD_rain"], self.settings["N_0rainD"], self.settings["n_moments"], self.settings["moments_file"],
		#input
		self.nlyr,self.ngridx,self.ngridy,self.nfreqs,self.freqs,
		self.year,self.month,self.day,self.time,
		self.deltax,self.deltay, self.lat,self.lon,self.model_i,self.model_j,
		self.wind10u,self.wind10v,self.lfrac,
		self.relhum_lev,self.press_lev,self.temp_lev,self.hgt_lev,
		self.iwv,self.cwp,self.iwp,self.rwp,self.swp,self.gwp,
		self.cwc_q,self.iwc_q,self.rwc_q,self.swc_q,self.gwc_q)
		
		self.r["Ze_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["attenuationHydro_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["attenuationAtmo_dimensions"] = ["gridx","gridy","lyr","frequency"]
		self.r["tb_dimensions"] = ["gridx","gridy","outlevels","angles","frequency","stokes"]
		
		self.r["settings"] = self.settings
		
		#for key in self.__dict__.keys():
			#print key
			#print self.__dict__[key]


def writeToNumpy(self,fname):
	try: 
		self.r
		self.r["pamtraVersion"]
	except:
		raise IOError ("run runPamtra first!")
		
	f = open(fname, "w")
	pickle.dump(self.r, f)
	f.close()

def loadFromNumpy(self,fname):
	try: 
		f = open(fname, "r")
		self.r = pickle.load(f)
		f.close()
	except:
		raise IOError ("Could not read data")


import pyPamtraLib
import numpy as np

class pyPamtra:

	
	def __init__(self,inputFile,frequency,userSettings):

		settings = dict()
		settings["verbose"]=0
		settings["write_nc"]=True
		settings["dump_to_file"]=False
		settings["tmp_path"]='/tmp/'
		settings["write_nc"]=False
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
		
		nlyr = 50
		ngridx = 4
		ngridy = 1
		freqs = [frequency,10]
		nfreqs = len(freqs)
		
		#output
		self.pamtraVersion,self.pamtraHash,\
		self.Ze,self.attenuationHydro,self.attenuationAtmo,self.hgt, self.angles = \
		pyPamtraLib.pypamtralib(inputFile,
		#settings
		settings["verbose"], settings["write_nc"], settings["dump_to_file"], settings["input_path"], settings["output_path"], settings["tmp_path"], settings["data_path"], settings["obs_height"], settings["units"], settings["outpol"], settings["freq_str"], settings["file_desc"], settings["creator"], settings["active"], settings["passive"], settings["ground_type"], settings["salinity"], settings["emissivity"], settings["lgas_extinction"], settings["gas_mod"], settings["lhyd_extinction"], settings["lphase_flag"], settings["SD_snow"], settings["N_0snowDsnow"], settings["EM_snow"], settings["SP"], settings["isnow_n0"], settings["liu_type"], settings["SD_grau"], settings["N_0grauDgrau"], settings["EM_grau"], settings["EM_ice"], settings["SD_rain"], settings["N_0rainD"], settings["n_moments"], settings["moments_file"],
		#input
		nlyr,ngridx,ngridy,nfreqs,freqs)
		
		for key in self.__dict__.keys():
			print key
			print self.__dict__[key]
		
		
		

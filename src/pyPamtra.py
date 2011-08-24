import pyPamtraLib

def pyPamtra(inputFile,frequency,userSettings):





	defaultSettings = dict()
	defaultSettings["verbose"]=0
	defaultSettings["write_nc"]=True
	defaultSettings["dump_to_file"]=False
	defaultSettings["tmp_path"]='/tmp/'
	defaultSettings["write_nc"]=False
	defaultSettings["input_path"]='../test/referenceProfile'
	defaultSettings["output_path"]='../test/tmp'
	defaultSettings["data_path"]='/home/mech/models/pamtra/data/'

	defaultSettings["obs_height"]=833000.
	defaultSettings["units"]='T'
	defaultSettings["outpol"]='VH'
	defaultSettings["freq_str"]=''
	defaultSettings["file_desc"]=''
	defaultSettings["creator"]='Pamtrauser'

	defaultSettings["active"]=True
	defaultSettings["passive"]=True

	defaultSettings["ground_type"]='S'
	defaultSettings["salinity"]=33.0
	defaultSettings["emissivity"]=0.6

	defaultSettings["lgas_extinction"]=True
	defaultSettings["gas_mod"]='R98'

	defaultSettings["lhyd_extinction"]=True
	defaultSettings["lphase_flag"]= True

	defaultSettings["SD_snow"]='Exp' 
	defaultSettings["N_0snowDsnow"]=7.628 
	defaultSettings["EM_snow"]='icesf' 
	defaultSettings["SP"]=0.2 
	defaultSettings["isnow_n0"]=1
	defaultSettings["liu_type"]=8

	defaultSettings["SD_grau"]='Exp' 
	defaultSettings["N_0grauDgrau"]=4.0 
	defaultSettings["EM_grau"]='surus'

	defaultSettings["EM_ice"]='mieic'

	defaultSettings["SD_rain"]='Exp' 
	defaultSettings["N_0rainD"]=8.0

	defaultSettings["n_moments"]=1
	defaultSettings["moments_file"]='snowCRYSTAL'

	for key in userSettings:
		try: defaultSettings[key] = userSettings[key]
		except: "Could not parse ",key
	
	pyPamtraLib.pypamtralib(inputFile,NamelistFile,frequency, defaultSettings["verbose"], defaultSettings["write_nc"], defaultSettings["dump_to_file"], defaultSettings["input_path"], defaultSettings["output_path"], defaultSettings["tmp_path"], defaultSettings["data_path"], defaultSettings["obs_height"], defaultSettings["units"], defaultSettings["outpol"], defaultSettings["freq_str"], defaultSettings["file_desc"], defaultSettings["creator"], defaultSettings["active"], defaultSettings["passive"], defaultSettings["ground_type"], defaultSettings["salinity"], defaultSettings["emissivity"], defaultSettings["lgas_extinction"], defaultSettings["gas_mod"], defaultSettings["lhyd_extinction"], defaultSettings["lphase_flag"], defaultSettings["SD_snow"], defaultSettings["N_0snowDsnow"], defaultSettings["EM_snow"], defaultSettings["SP"], defaultSettings["isnow_n0"], defaultSettings["liu_type"], defaultSettings["SD_grau"], defaultSettings["N_0grauDgrau"], defaultSettings["EM_grau"], defaultSettings["EM_ice"], defaultSettings["SD_rain"], defaultSettings["N_0rainD"], defaultSettings["n_moments"], defaultSettings["moments_file"])


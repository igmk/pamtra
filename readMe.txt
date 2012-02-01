Pamtra and pyPamtra have the following options.

Pamtra uses the nml file, pyPamtra the pyPamtra.set dictionary

verbose_mode:
	verbose: level of verbosity (int 0-?, default 0)

inoutput_mode:
	write_nc: write results to netcdf file instead of ASCII, pamtra only (bool, default true)
	dump_to_file: write old input files for RT3 (for debugging only, bool, default false)
	input_path: pamtra only (str, default profiles)
	output_path: pamtra only (str, default output)
	data_path: path containing the surface reflectivity data etc. (str, default data)
	tmp_path: (default tmp)

output
	obs_height=833000. 
	units='T'
	outpol='VH'
	freq_str=''
	file_desc=''
	creator: for netcdf file (str, default "Pamtra user")
	zeSplitUp: save Ze and Attenuation seperately for every hydrometeor kind, pyPamtra outputs both, option only important for writeResultsToNetcdf! (bool, default true)

run_mode
	active: calculate Ze and Attenuation (bool, default true)
	passive: calculate TB, thus run RT3 (bool, default true)

surface_params
	ground_type='S'
	salinity=33.0
	emissivity=0.6

gas_abs_mod
	lgas_extinction=.true.
	gas_mod='R98'

hyd_opts
	lhyd_extinction=.true.
	lphase_flag = .true.

snow_params
	SD_snow='Exp' 
	N_0snowDsnow=7.628 
	EM_snow='surus' 
	SP=0.2 
	isnow_n0=1
	liu_type=8

graupel_params
	SD_grau='Exp' 
	N_0grauDgrau=4.0 
	EM_grau='surus'

ice_params
	EM_ice='icesf'

rain_params
	SD_rain='Exp' 
	N_0rainD=8.0

moments
	n_moments=1
	moments_file='snowCRYSTAL'

pyPamtra only:
	freqs: used frequencies for calculations (list, default empty, in Pamtra handled by command line parameter)
	nfreqs: amount of frequencies (int, default 0)
	pyVerbose: Verbosity of the python part (int, default 0)
	
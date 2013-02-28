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

radar_simulator
	
	radar_nfft: number of FFT points in the Doppler spectrum [typically 256 or 512] (default 256)
	radar_no_Ave_ number of average spectra for noise variance reduction, typical range [1 150] (default 150)
	radar_max_V:MinimumNyquistVelocity in m/sec (default 7.885)
	radar_min_V:MaximumNyquistVelocity in m/sec (default -7.885)
	radar_turbulence_st: turbulence broadening standard deviation st, typical range [0.1 - 0.4] m/sec (default 0.15)
	radar_pnoise: radar noise in same unit as Ze mm⁶/m³ (default 1.d-3)

	radar_airmotion: is teh air in the radar volume moving vertically? (default  .false.)
	radar_airmotion_model: constant air movmement or non uniform beam filling: linear or step function ["constant","linear","step"] (default  "step")
	radar_airmotion_vmin: for nun uniform beamfilling minimal velocity, also taken for constant air movment(default  -4.d0)
	radar_airmotion_vmax: for nun uniform beamfilling minimal velocity, ignored for constant air movment(default  +4.d0)
	radar_airmotion_linear_steps: no of steps for linear approximation(default 30)
	radar_airmotion_step_vmin: ratio of volume which moves with vmin for step function(default  0.5d0)
	
	models available for fall velocity approximation
	radar_fallVel_cloud: (default "khvorostyanov01_drops")
	radar_fallVel_rain: (default  "khvorostyanov01_drops")
	radar_fallVel_ice: (default "khvorostyanov01_particles")
	radar_fallVel_snow: (default "khvorostyanov01_particles")
	radar_fallVel_graupel: (default "khvorostyanov01_spheres")
	radar_fallVel_hail: (default "khvorostyanov01_spheres")

	radar_aliasing_nyquist_interv: simulate aliasing effects: spectrum is added x time to the left and right.(default  1)
	radar_save_noise_corrected_spectra: for debugging: save the radar spectrum with noise removed (default  .false.)
	radar_use_hildebrand: use hildebrand & sekhon for noise estimation, actually not needed since noise is added artifically (default  .false.)
	radar_min_spectral_snr: threshold for peak detection. if radar_no_Ave >> 150, it can be set to 1.1(default  1.2)
	radar_convolution_fft: use fft for convolution of spectrum. is alomst 10 times faster, but can introduce aretfacts for radars with *extremely* low noise levels or if noise is turned off at all.  (default  .true.)


pyPamtra only:
	freqs: used frequencies for calculations (list, default empty, in Pamtra handled by command line parameter)
	nfreqs: amount of frequencies (int, default 0)
	pyVerbose: Verbosity of the python part (int, default 0)
	
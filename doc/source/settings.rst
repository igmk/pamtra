
Settings/Nml-File
=================

This file defines general settings of the PAMTRA run (e.g. whether active or passive simulations, just radar moments or full spectrum, paths for in- and output, etc.). The full information of available options can be found in settings.f90

================================== ============================== =========== ===========================================================================================
Variable                           Values                         Default     Description                                                                               
================================== ============================== =========== ===========================================================================================
active                             bool                           True        Activate radar simulator
creator                            str                            Pamtrauser  Netcdf file creator
data_path                          str                            data/              
emissivity                         positive float                 0.6
file_desc                          str                            ""
gas_mod                            str                            R98
ground_type                        str                            L
hydro_adaptive_grid                bool                           True
hydro_fullspec                     bool                           False
hydro_includehydroinrhoair         bool                           True
hydro_limit_density_area           bool                           True
hydro_softsphere_min_density       positive float                 10.0
hydro_threshold                    positive float                 1e-10       minimum required hydrometeor concentration kg/m^3
lgas_extinction                    bool                           True        gas extinction desired
lhyd_extinction                    bool                           True        hydrometeor extinction desired
liq_mod                            str                            Ell
obs_height                         positive float                 833000.0    upper level output height [m] (> 100000. for satellite)
outpol                             str                            VH
passive                            bool                           True        estimate brightness temperatures
radar\_airmotion                   boolean                        False       Consider air motion in direction of radar beam.                                                                                                                                                                                                                        
radar\_airmotion\_linear\_steps    positive integer               30          For linear function: number of discrete intervals.                                                                                                                                                                                                                     
radar\_airmotion\_model            constant, linear, step         step        Model to describe vertical air motion: Either constant velocity, linear change from vmin to vmax or abrupt change using a step function.                                                                                                                               
radar\_airmotion\_step\_vmin       positive float                 0.5         For step function: volume ratio between vmin and vmax.                                                                                                                                                                                                                 
radar\_airmotion\_vmin             float                          -4 m/s      Minimal air motion of for step and linear function. Also used for constant air motion.                                                                                                                                                                                 
radar\_airmotion\_vmax             float                          4 m/s       Maximal air motion of for step and linear function.                                                                                                                                                                                                                    
radar_aliasing_nyquist_interv      positive integer               1
radar\_attenuation                 disabled, bottom-up, top-down  disabled    Attenuate radar spectrum and  Z_e  depending on measurement geometry
radar_convolution_fft              boolean                        True
radar\_K2 (|K_w^2|)                positive float                 0.93        Dielectric factor of water
radar\_max\_v ( v_nyq )            float                          -7.885 m/s  Maximum Nyquist velocity (usually radar\_min\_V = -radar\_max\_V)   
radar\_min\_v ( v_nyq )            float                          7.885 m/s   Minimum Nyquist velocity 
radar_min_spectral_snr             positive float                 1.2
radar\_mode                        simple, spectrum, moments      simple      Use "simple" radar simulator provides only Z_e by integrating Eq. \label{eq:etaD} over  D. The advanced "spectrum" simulator simulates the complete radar Doppler spectrum and estimates all moments from the spectrum. "moments" is identical to "spectrum" but the full Doppler spectrum is discarded to save memory. 
radar\_nfft ( N_fft )              positive integer               256         Number of FFT points in the Doppler spectrum 
radar\_no\_Ave ( Nave )            positive integer               150         Number of spectral averages
radar_noise_distance_factor        positive float                 2.0
radar_npeaks                       positive integer               1
radar\_pnoise0 ( N_1000 )          float                          -32.23 dBz  Radar noise at 1km in same unit as reflectivity Z_e (Eq.~\ref{eq:radarnoise})
radar\_polarisation                NN, HV, VH, VV, HH             NN          Radar polarisation. NN: no polarisation, HV: horizontal transmit, vertical receive, etc.. Can be a comma separated list.
radar_receiver_uncertainty_std     positive float                 0.0
radar_save_noise_corrected_spectra boolean                        False
radar_smooth_spectrum              boolean                        True        smooth spectrum before estimating moments
radar\_use\_hildebrand             boolean                        False       Derive  N_P  not from radar\_pnoise0 but using the method of \citet{hildebrand:1974a}.                                                                                                          
randomseed                         integer                        0           0 is real noise, -1 means that the seed is created from latitude and longitude, other value gives always the same random numbers
salinity                           float                          33.0        sea surface salinity
save_psd                           boolean                        False       also saves the PSDs used for radiative transfer
save_ssp                           boolean                        False       also saves the single scattering properties used for radiative transfer
tmatrix_db                         none or file                   none        use data base to cache T-Matrix calculations
tmatrix_db_path                    str                            database/   path to T-Matrix data base
write_nc                           bool                           True        write netcdf or ascii output
================================== ============================== =========== ===========================================================================================



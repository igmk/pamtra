..  _settings:


Settings/Nml-File
=================


This file defines general settings of the PAMTRA run (e.g., whether active or passive simulations, just radar moments or full spectrum, paths for in- and output, etc.). These settings are stored in a dictionary *nmlSet* of the pyPamtra object. In addition, some low level settings, like verbosity of FORTRAN and python, are stored in the *set* dictionary. The full information of available options can be found in *settings.f90* in the PAMTRA source directory.

nmlSet() settings
*****************
================================== ============================== ================== =============================================================================================================================================================================================================================================================================================================================================================================
Variable                           Values                         Default            Description
================================== ============================== ================== =============================================================================================================================================================================================================================================================================================================================================================================
active                             bool                           True               Activate radar simulator
add_obs_height_to_layer            bool                           False              If observation heights for the output are selected, these can be added as additional atmospheric layer boudnaries. In case the vertical grid of observation layers is very fine, it can happen that several heights are within the same layer and therefore give the same result. For cases where the observation height is very far away from the next atmospheric model layer, it might not be 100% representative and it could be beneficial as well.
conserve_mass_rescale_dsd          bool                           True               In case the mass mixing ratio for an hydrometeor calculated integrating the drop-size-distribution (DSD) doesn't correspond to the input value, rescale the DSD to account for the mass loss.
creator                            str                            Pamtrauser         Netcdf file creator
data_path                          str                            $PAMTRA_DATADIR    Path for emissivity files and other data. If value is $PAMTRA_DATADIR, the corresponding environment variable is used.
emissivity                         positive float [0,1]           0.6                Surface emissivity used for both polarizations
file_desc                          str                            ""                 In pure FORTRAN mode and netCDF output, this string is used as an extension to the output file name. For sensitivity studies this might be helpful. 
gas_mod                            L93, R98                       R98                Model for gas absorption. Either ROSENKRANZ (R98) or LIEBE (L93)
hydro_adaptive_grid                bool                           True
hydro_fullspec                     bool                           False              For pyPamtra only: Do not estimate particle diameter, mass, area, number concentration, rho and aspect ratio directly from the descriptor file but pass them directly from python to PAMTRA using numpy arrays. See also addFullSpectra() of pyPamtra's descriptorFile class.
hydro_includehydroinrhoair         bool                           True               Include hydrometeors when estimating the density of wet air. Different models use different conventions here.
hydro_limit_density_area           bool                           True               Change mass, cross section area and density of particles in case it is larger or smaller than possible. Min density is hydro_softsphere_min_density, max density is 917 kg/m\ :sup:`3`. max area is D\ :sup:`2`
hydro_softsphere_min_density       positive float                 10.0               If hydro_limit_density_area=True, limit minimal density to this value.
hydro_threshold                    positive float                 1e-10              minimum required hydrometeor concentration kg/m\ :sup:`3`.
lgas_extinction                    bool                           True               gas extinction desired
lhyd_extinction                    bool                           True               hydrometeor extinction desired
liq_mod                            str                            Ell
obs_height                         positive float                 833000.0           upper level output height [m] (> 100000. for satellite)
outpol                             str                            VH
passive                            bool                           True               estimate brightness temperatures
radar_allow_negative_dD_dU         bool                           False              allow that particle velocity is decreasing with size. Should be usually set to false.
radar\_airmotion                   boolean                        False              Consider air motion in direction of radar beam.
radar\_airmotion\_linear\_steps    positive integer               30                 For linear function: number of discrete intervals.
radar\_airmotion\_model            constant, linear, step         step               Model to describe vertical air motion: Either constant velocity, linear change from vmin to vmax or abrupt change using a step function.
radar\_airmotion\_step\_vmin       positive float                 0.5                For step function: volume ratio between vmin and vmax.
radar\_airmotion\_vmin             float                          -4 m/s             Minimal air motion of for step and linear function. Also used for constant air motion.
radar\_airmotion\_vmax             float                          4 m/s              Maximal air motion of for step and linear function.
radar_aliasing_nyquist_interv      positive integer               1                  Consider aliasing effects for overspending the nyquist range radar_aliasing_nyquist_interv times.
radar\_attenuation                 disabled, bottom-up, top-down  disabled           Attenuate radar spectrum and  Z_e  depending on measurement geometry (bottom-up for upward looking, top-down for downward-looking).
radar_convolution_fft              boolean                        True               Use FFT for convolution. FFt is much faster, but can have numerical issues in rare cases.
radar_fwhr_beamwidth_deg           float*                         0.3                radar full width half radiation beamwidth (required for spectral broadening estimation)
radar_integration_time             float*                         1.4                radar beamwidth (required for spectral broadening estimation)
radar\_K2 (\|K_w^2\|)                positive float*                0.93               Dielectric factor of water used to estimate radr reflectivity.
radar\_max\_v ( v_nyq )            float*                         -7.885 m/s         Maximum Nyquist velocity (usually radar\_min\_V = -radar\_max\_V)
radar\_min\_v ( v_nyq )            float*                         7.885 m/s          Minimum Nyquist velocity
radar_peak_min_bins                int*                           2                  Minimum peak width
radar_peak_min_snr                 float*                         -10 dB             Minimal required SNR reqired for a peak. See radar_peak_min_snr for defintion
radar_peak_snr_definition          specLin | log                  log                log: radar_peak_min_snr describes snr of peak in dB. linSpec: radar_peak_min_snr descibes mean signal+noise to noise ratio (available for historical reasons)
radar\_mode                        simple, spectrum, moments      simple             Use "simple" radar simulator provides only Z_e by integrating over D. The advanced "spectrum" simulator simulates the complete radar Doppler spectrum and estimates all moments from the spectrum. "moments" is identical to "spectrum" but the full Doppler spectrum is discarded to save memory.
radar\_nfft ( N_fft )              positive integer               256                Number of FFT points in the Doppler spectrum
radar\_no\_Ave ( Nave )            positive integer*              150                Number of spectral averages
radar_noise_distance_factor        positive float*                2.0                Required distance of the peak edge to the noise level. If radar_noise_distance_factor<0 and radar\_use\_hildebrand, then noise_max from Hildebrand is used for peak edge determination. Sometimes, lower SNR values can be achieved with radar_noise_distance_factor instead of noise_max
radar_npeaks                       1                              1                  Number of detected peaks in the Doppler spectrum. As of today fixed to 1.
radar\_pnoise0 ( N_1000 )          float*                         -32.23 dBz         Radar noise at 1km in same unit as reflectivity Z_e
radar\_polarisation                NN, HV, VH, VV, HH             NN                 Radar polarisation. NN: no polarisation, HV: horizontal transmit, vertical receive, etc.. Can be a comma separated list.
radar_receiver_miscalibration      float*                         0.0 dB             Radar calibration error
radar_receiver_uncertainty_std     positive float*                0.0                Add Gaussian noise to radar noise level to simulate unstable receivers
radar_save_noise_corrected_spectra boolean                        False              For debugging purposes: Save radar Doppler spectrum after noise is removed
radar_smooth_spectrum              boolean                        True               smooth spectrum before estimating moments
radar\_use\_hildebrand             boolean                        False              Derive  N_P  not from radar\_pnoise0 but using the method of \citet{hildebrand:1974a}. Set  radar_noise_distance_factor<0 to use also noise_max from hildebrand for determination od the peak edge. Sometimes, lower SNR values can be achieved with radar_noise_distance_factor instead of noise_max
radar_use_wider_peak               boolean                        False              Include the found peak edge (if peak edge is still larger than mean noise) into the peak which is used for moment estimation.
randomseed                         integer                        0                  0 is real noise, -1 means that the seed is created from latitude and longitude, other value gives always the same random numbers
read_turbulence_ascii              bool                           False              If .true. turbulence need to be included in the ascii input_file, rightmost column. Not relevant for pyPamtra and for passive simulations.
salinity                           float                          33.0               sea surface salinity
save_psd                           boolean                        False              also saves the PSDs used for radiative transfer
save_ssp                           boolean                        False              also saves the single scattering properties used for radiative transfer
tmatrix_db                         none or file                   none               use data base to cache T-Matrix calculations
tmatrix_db_path                    str                            database/          path to T-Matrix data base
write_nc                           bool                           True               write netcdf or ascii output
================================== ============================== ================== =============================================================================================================================================================================================================================================================================================================================================================================

\* These variables *can* be also provided as list to account for different instrument specifications. In this case, each entry corresponds to one frequency.

set() settings
**************
================== ============================== ================== ==============================================================================================================================================================================================================================================================================================================================================================================
Variable           Values                         Default            Description
================== ============================== ================== ==============================================================================================================================================================================================================================================================================================================================================================================
verbose            positive integer               0                  Verbosity of the FORTRAN routines
pyVerbose          positive integer               0                  Verbosity of the pyPamtra python modules
namelist_file      str                            TMPFILE            path and name of the FORTRAN namelist file
freqs              list of float                  empty              list of frequencies, set automatically at program start
================== ============================== ================== ==============================================================================================================================================================================================================================================================================================================================================================================

Other default settings

================== ============================== ================== ==============================================================================================================================================================================================================================================================================================================================================================================
Variable           Values                         Default            Description
================== ============================== ================== ==============================================================================================================================================================================================================================================================================================================================================================================
sfc_refl           S,L,F                          S                  Specular, Lambertian, or Fresnel
================== ============================== ================== ==============================================================================================================================================================================================================================================================================================================================================================================

module nml_params
  ! Description:
  ! Definition of all name list paramters for pamtra
  !and global settings!


  use kinds

  implicit none

  !!Global Stettings
  integer, parameter :: MAXV = 64,   &
       MAXLAY = 600, &
       maxleg = 200, &
       maxfreq = 100, &
       nummu = 16, & ! no. of observation angles
       NSTOKES = 2, &
       NOUTLEVELS = 2
  integer, parameter :: SRC_CODE = 2,&
       NUMAZIMUTHS = 1,&
       Aziorder = 0

  real(kind=dbl), parameter :: SKY_TEMP    = 2.73d0,  &    ! cosmic background
       DIRECT_FLUX = 0.d0 ,&
       DIRECT_MU   = 0.0d0 

  character(1), parameter :: QUAD_TYPE = 'L',&
       DELTAM = 'N'

  !!Set by namelist file
  integer :: verbose, n_moments, isnow_n0, liu_type

  real(kind=dbl) :: obs_height     ! upper level output height [m] (> 100000. for satellite)
  real(kind=dbl) :: emissivity
  real(kind=dbl) :: N_0rainD, N_0snowDsnow, N_0grauDgrau, N_0hailDhail,  SP
  real(kind=dbl) :: as_ratio, snow_density, graupel_density, hail_density
  real(kind=dbl) :: salinity         ! sea surface salinity
  integer :: radar_nfft !number of FFT points in the Doppler spectrum [typically 256 or 512]
  integer :: radar_no_Ave !number of average spectra for noise variance reduction, typical range [1 40]
  integer :: radar_airmotion_linear_steps !for linear air velocity model, how many staps shall be calculated?
  integer :: radar_aliasing_nyquist_interv !how many additional nyquists intervalls shall be added to the spectrum to deal with aliasing effects
  real(kind=dbl) :: radar_max_V !MinimumNyquistVelocity in m/sec
  real(kind=dbl) :: radar_min_V !MaximumNyquistVelocity in m/sec
  real(kind=dbl) :: radar_turbulence_st !turbulence broadening standard deviation st, typical range [0.1 - 0.4] m/sec
  real(kind=dbl) :: radar_pnoise !radar noise
  real(kind=dbl) :: radar_airmotion_vmin
  real(kind=dbl) :: radar_airmotion_vmax
  real(kind=dbl) :: radar_airmotion_step_vmin

  logical ::  in_python !are we in python

  logical :: dump_to_file, &   ! flag for profile and ssp dump
       lphase_flag, &        ! flag for phase function calculation
       lgas_extinction, &    ! gas extinction desired
       lhyd_extinction, &    ! hydrometeor extinction desired
       use_rain_db, &    ! use the tmatrix database for rain
       use_snow_db, &    ! use the tmatrix database for snow
       write_nc, &  ! write netcdf or ascii output
       active, &  	   ! calculate active stuff
       passive, &     ! calculate passive stuff (with RT3)
       activeLogScale, & !save ze and att in log scale or linear
       jacobian_mode, &  ! special jacobian mode which does not calculate the whole scattering properties each time. only rt4!
       radar_airmotion, &   ! apply vertical air motion
       radar_save_noise_corrected_spectra, & !remove the noise from the calculated spectrum again (for testing) 
       radar_use_hildebrand  ! use Hildebrand & Sekhon for noise estimation as a real radar would do. However, since we set the noise (radar_pnoise) we can skip that.

  character(5) :: EM_ice, EM_snow, EM_grau, EM_hail, EM_cloud, EM_rain
  character(1) :: SD_cloud, SD_ice, SD_rain, SD_snow, SD_grau, SD_hail
  character(3) :: gas_mod
  character(20) :: moments_file,file_desc
  character(100) :: input_path, output_path, tmp_path,creator, data_path
  character(18) :: freq_str
  character(2) :: OUTPOL
  character(1) :: GROUND_TYPE, UNITS
  character(3) :: rt_mode
  character(10) :: input_type, crm_case
  character(100) :: crm_data, crm_data2, crm_constants
  character(8) :: radar_airmotion_model, radar_mode
  character(30) :: radar_fallVel_cloud, radar_fallVel_rain, radar_fallVel_ice,&
		    radar_fallVel_snow, radar_fallVel_graupel, radar_fallVel_hail

  integer :: radar_nfft_aliased, radar_maxTurbTerms !are gained from radar_aliasing_nyquist_interv and radar_nfft
contains

  subroutine nml_params_read

    use file_mod, only: namelist_file

    ! name list declarations
    namelist / verbose_mode / verbose
    namelist / inoutput_mode / input_path, output_path,&
         tmp_path, dump_to_file, write_nc, data_path,&
         input_type, crm_case, crm_data, crm_data2, crm_constants, &
	 jacobian_mode
    namelist / output / obs_height,units,outpol,freq_str,file_desc,creator, &
	  activeLogScale
    namelist / run_mode / active, passive,rt_mode, radar_mode
    namelist / surface_params / ground_type,salinity, emissivity
    namelist / gas_abs_mod / lgas_extinction, gas_mod
    namelist / hyd_opts / lhyd_extinction, lphase_flag
    namelist / cloud_params / SD_cloud, EM_cloud
    namelist / ice_params / SD_ice, EM_ice
    namelist / rain_params / SD_rain, N_0rainD, use_rain_db, EM_rain
    namelist / snow_params / SD_snow, N_0snowDsnow, EM_snow, use_snow_db, as_ratio,snow_density, SP, isnow_n0, liu_type
    namelist / graupel_params / SD_grau, N_0grauDgrau, EM_grau, graupel_density
    namelist / hail_params / SD_hail, N_0hailDhail, EM_hail, hail_density
    namelist / moments / n_moments, moments_file
    namelist / radar_simulator / radar_nfft,radar_no_Ave, radar_max_V, radar_min_V, &
	       radar_turbulence_st, radar_pnoise, radar_airmotion, radar_airmotion_model, &
	       radar_airmotion_vmin, radar_airmotion_vmax, radar_airmotion_linear_steps, &
	       radar_airmotion_step_vmin, radar_fallVel_cloud, radar_fallVel_rain, radar_fallVel_ice,&
	       radar_fallVel_snow, radar_fallVel_graupel, radar_fallVel_hail, radar_aliasing_nyquist_interv, &
	       radar_save_noise_corrected_spectra, radar_use_hildebrand

    !set namelist defaults!
    ! sec verbose_mode
    verbose=0
    ! sec inoutput_mode
    write_nc=.true.
    dump_to_file=.false.
    input_path='profile/'
    output_path='output/'
    input_type='profile'
    tmp_path='/tmp/'
    data_path='data/'
    crm_case=''
    crm_data=''
    crm_data2=''
    crm_constants=''
    jacobian_mode=.false. !profile 1,1 is reference, for all other colums only layers with different values are calculated
    ! sec output
    obs_height=833000.
    units='T'
    outpol='VH'
    freq_str=''
    file_desc=''
    creator='Pamtrauser'
    activeLogScale = .true.
    ! sec run_mode
    active=.true.
    passive=.true.
    rt_mode='rt3'
    radar_mode="simple" !"splitted"|"moments"|"spectrum"
    ! sec surface params
    ground_type='S'
    salinity=33.0
    emissivity=0.6
    ! sec gas_abs_mod
    lgas_extinction=.true.
    gas_mod='R98'
    ! sec hyd_opts
    lhyd_extinction=.true.
    lphase_flag = .true.
    ! sec cloud_params
    SD_cloud='C'
    EM_cloud="miecl"
    ! sec ice_params
    SD_ice='C'
    EM_ice='mieic'
    ! sec rain_params
    SD_rain='C'
    N_0rainD=8.0
    use_rain_db=.true.
    EM_rain="miera"
    ! sec snow_params
    SD_snow='C'
    N_0snowDsnow=7.628 
    EM_snow='densi'
    use_snow_db=.true.
    as_ratio=0.5d0
    snow_density=200.d0
    SP=0.2 
    isnow_n0=1
    liu_type=8
    ! sec graupel_params
    SD_grau='C'
    N_0grauDgrau=4.0 
    EM_grau='densi'
    graupel_density=400.d0
    ! sec hail_params
    SD_hail='C'
    N_0hailDhail=4.0
    EM_hail='densi'
    hail_density=917.d0
    ! sec moments
    n_moments=1
    moments_file='snowCRYSTAL'
    ! radar_simulator
    !number of FFT points in the Doppler spectrum [typically 256 or 512]
    radar_nfft=256
    !number of average spectra for noise variance reduction, typical range [1 150]
    radar_no_Ave=150
    !MinimumNyquistVelocity in m/sec
    radar_max_V=7.885
    !MaximumNyquistVelocity in m/sec
    radar_min_V=-7.885
    !turbulence broadening standard deviation st, typical range [0.1 - 0.4] m/sec
    radar_turbulence_st=0.15
      !radar noise in same unit as Ze mm⁶/m³
    radar_pnoise=1.d-3

    radar_airmotion = .true.
    radar_airmotion_model = "step" !"constant","linear","step"
    radar_airmotion_vmin = -4.d0
    radar_airmotion_vmax = +4.d0
    radar_airmotion_linear_steps = 30
    radar_airmotion_step_vmin = 0.5d0

    radar_fallVel_cloud ="khvorostyanov01_drops"
    radar_fallVel_rain = "khvorostyanov01_drops"
    radar_fallVel_ice ="khvorostyanov01_particles"
    radar_fallVel_snow ="khvorostyanov01_particles"
    radar_fallVel_graupel ="khvorostyanov01_spheres"
    radar_fallVel_hail ="khvorostyanov01_spheres"

    radar_aliasing_nyquist_interv = 1
    radar_save_noise_corrected_spectra = .false.
    radar_use_hildebrand = .false.

    ! read name list parameter file
    open(7, file=namelist_file,delim='APOSTROPHE')
    read(7,nml=verbose_mode)
    read(7,nml=inoutput_mode)
!    if (verbose .gt. 1) print*, input_path, output_path,tmp_path, dump_to_file, write_nc, data_path
    read(7,nml=output)
!    if (verbose .gt. 1) print*, obs_height,units,outpol,freq_str,file_desc,creator,zeSplitUp
    read(7,nml=run_mode)
!    if (verbose .gt. 1) print*, active, passive,rt_mode
    read(7,nml=surface_params)
!    if (verbose .gt. 1) print*, ground_type,salinity, emissivity
    read(7,nml=gas_abs_mod)
!    if (verbose .gt. 1) print*, lgas_extinction, gas_mod
    read(7,nml=hyd_opts)
!    if (verbose .gt. 1) print*, lhyd_extinction, lphase_flag
    read(7,nml=cloud_params)
!    if (verbose .gt. 1) print*, SD_cloud
    read(7,nml=ice_params)
!    if (verbose .gt. 1) print*, SD_ice, EM_ice
    read(7,nml=rain_params)
!    if (verbose .gt. 1) print*, SD_rain, N_0rainD, use_rain_db
    read(7,nml=snow_params)
!    if (verbose .gt. 1) print*, SD_snow, N_0snowDsnow, EM_snow, use_snow_db, as_ratio,snow_density, SP, isnow_n0, liu_type
    read(7,nml=graupel_params)
!    if (verbose .gt. 1) print*, SD_grau, N_0grauDgrau, EM_grau, graupel_density
    read(7,nml=hail_params)
!    if (verbose .gt. 1) print*, SD_hail, N_0hailDhail, EM_hail, hail_density
    read(7,nml=moments)
!    if (verbose .gt. 1) print*, n_moments, moments_file
    read(7,nml=radar_simulator)

    close(7)

!mix some variables to make new ones:
  radar_nfft_aliased = radar_nfft *(1+2*radar_aliasing_nyquist_interv)
  radar_maxTurbTerms = radar_nfft_aliased * 12

     if (verbose .gt. 3) then  
    print*, "verbose_mode ",  verbose
    print*, "inoutput_mode ",  input_path, output_path,&
         tmp_path, dump_to_file, write_nc, data_path,&
         input_type, crm_case, crm_data, crm_data2, crm_constants, &
	 jacobian_mode
    print*, "output ",  obs_height,units,outpol,freq_str,file_desc,creator, &
	  activeLogScale
    print*, "run_mode ",  active, passive,rt_mode, radar_mode
    print*, "surface_params ",  ground_type,salinity, emissivity
    print*, "gas_abs_mod ",  lgas_extinction, gas_mod
    print*, "hyd_opts ",  lhyd_extinction, lphase_flag
    print*, "cloud_params ",  SD_cloud, EM_cloud
    print*, "ice_params ",  SD_ice, EM_ice
    print*, "rain_params ",  SD_rain, N_0rainD, use_rain_db, EM_rain
    print*, "snow_params ",  SD_snow, N_0snowDsnow, EM_snow, use_snow_db, as_ratio,snow_density, SP, isnow_n0, liu_type
    print*, "graupel_params ",  SD_grau, N_0grauDgrau, EM_grau, graupel_density
    print*, "hail_params ",  SD_hail, N_0hailDhail, EM_hail, hail_density
    print*, "moments ",  n_moments, moments_file
    print*, "radar_simulator ",  radar_nfft,radar_no_Ave, radar_max_V, radar_min_V, &
	       radar_turbulence_st, radar_pnoise, radar_airmotion, radar_airmotion_model, &
	       radar_airmotion_vmin, radar_airmotion_vmax, radar_airmotion_linear_steps, &
	       radar_airmotion_step_vmin, radar_save_noise_corrected_spectra, radar_use_hildebrand

  end if

    if (n_moments .ne. 1 .and. n_moments .ne. 2) stop "n_moments is not 1 or 2"

    if (verbose .gt. 1) print *,"PASSIVE: ", passive, "ACTIVE: ", active


    return
  end subroutine nml_params_read
end module nml_params

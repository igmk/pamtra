module settings
    ! Description:
    ! Definition of all command line and name list paramters for pamtra
    ! and global settings!

    use kinds
    use report_module
    implicit none
    save

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

    real(kind=dbl), parameter :: DIRECT_FLUX = 0._dbl ,&
    DIRECT_MU   = 0._dbl

    character(1), parameter :: QUAD_TYPE = 'L',&
    DELTAM = 'N'

    ! set by command line options
    integer(kind=long) :: nfrq

    !!Set by namelist file
    integer(kind=long):: n_moments

    real(kind=dbl) :: obs_height     ! upper level output height [m] (> 100000. for satellite)
    real(kind=dbl) :: emissivity
    real(kind=dbl) :: salinity         ! sea surface salinity
!     double precision, dimension(maxfreq) :: freqs
    real(kind=dbl), dimension(maxfreq) :: freqs

    integer(kind=long) :: radar_nfft !number of FFT points in the Doppler spectrum [typically 256 or 512]
    integer(kind=long):: radar_no_Ave !number of average spectra for noise variance reduction, typical range [1 40]
    integer(kind=long):: radar_airmotion_linear_steps !for linear air velocity model, how many staps shall be calculated?
    integer(kind=long):: radar_aliasing_nyquist_interv !how many additional nyquists intervalls shall be added to the spectrum to deal with aliasing effects
    integer(kind=long) :: radar_nPeaks
    real(kind=dbl) :: radar_max_V !MinimumNyquistVelocity in m/sec
    real(kind=dbl) :: radar_min_V !MaximumNyquistVelocity in m/sec
    real(kind=dbl) :: radar_pnoise0 !radar noise at 1km
    real(kind=dbl) :: radar_airmotion_vmin
    real(kind=dbl) :: radar_airmotion_vmax
    real(kind=dbl) :: radar_airmotion_step_vmin
    real(kind=dbl) :: radar_min_spectral_snr !threshold for peak detection
    real(kind=dbl) :: radar_K2
    real(kind=dbl) :: radar_receiver_uncertainty_std
    real(kind=dbl) :: hydro_softsphere_min_density !tmatrix method numerically unstable for extremely low density

    real(kind=dbl) :: hydro_threshold, radar_noise_distance_factor

  integer, parameter :: maxnleg = 200 !max legnth of legendre series
  
    logical :: in_python !are we in python

    logical :: dump_to_file, &   ! flag for profile and ssp dump
    lphase_flag, &        ! flag for phase function calculation
    lgas_extinction, &    ! gas extinction desired
    lhyd_extinction, &    ! hydrometeor extinction desired
    write_nc, &  ! write netcdf or ascii output
    active, &  	   ! calculate active stuff
    passive, &     ! calculate passive stuff (with RT4)
    jacobian_mode, &  ! special jacobian mode which does not calculate the whole scattering properties each time. only rt4!
    radar_airmotion, &   ! apply vertical air motion
    radar_save_noise_corrected_spectra, & !remove the noise from the calculated spectrum again (for testing)
    radar_use_hildebrand,&  ! use Hildebrand & Sekhon for noise estimation as a real radar would do. However, since we set the noise (radar_pnoise0) we can skip that.
    radar_convolution_fft,&!use fft for convolution of spectrum
    save_psd, &
    radar_smooth_spectrum, &
    hydro_fullSpec, &
    hydro_limit_density_area, &
    add_obs_height_to_layer

    character(3) :: gas_mod
    character(20) :: moments_file,file_desc
    character(100) :: input_path, output_path, tmp_path,creator, data_path
    character(18) :: freq_str
    character(2) :: OUTPOL
    character(1) :: GROUND_TYPE, UNITS
    character(10) :: input_type, crm_case
    character(100) :: crm_data, crm_data2, crm_constants
    character(8) :: radar_airmotion_model, radar_mode

    character(99)  :: input_file        ! name of profile
    character(300) :: namelist_file     ! name of nml_file
    character(300) :: nc_out_file       ! name of netcdf output file
    character(9) :: frq_str_s,frq_str_e
    character(8), dimension(maxfreq) :: frqs_str
    character(300) :: descriptor_file_name

    integer(kind=long):: radar_nfft_aliased, radar_maxTurbTerms !are gained from radar_aliasing_nyquist_interv and radar_nfft
    
    integer(kind=long) :: randomseed !random seed, 0 means time dependence
contains

    subroutine settings_read(errorstatus)

    use kinds
    implicit none
    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'settings_read'

        ! name list declarations
        namelist / inoutput_mode / input_path, output_path,&
        tmp_path, dump_to_file, write_nc, data_path,&
        input_type, crm_case, crm_data, crm_data2, crm_constants, &
        jacobian_mode, save_psd
        namelist / output / obs_height,units,outpol,freq_str,file_desc,creator, add_obs_height_to_layer
        namelist / run_mode / active, passive,radar_mode, randomseed
        namelist / surface_params / ground_type,salinity, emissivity
        namelist / gas_abs_mod / lgas_extinction, gas_mod
        namelist / hyd_opts / lhyd_extinction, lphase_flag, hydro_fullSpec, hydro_limit_density_area,&
                  hydro_softsphere_min_density
	namelist / moments / n_moments, moments_file
	namelist / radar_simulator / radar_nfft,radar_no_Ave, radar_max_V, radar_min_V, &
		  radar_pnoise0, radar_airmotion, radar_airmotion_model, &
		  radar_airmotion_vmin, radar_airmotion_vmax, radar_airmotion_linear_steps, &
		  radar_airmotion_step_vmin, radar_aliasing_nyquist_interv, &
		  radar_save_noise_corrected_spectra, radar_use_hildebrand, radar_min_spectral_snr, radar_convolution_fft, &
                  radar_K2, radar_noise_distance_factor, radar_receiver_uncertainty_std,&
                  radar_nPeaks, radar_smooth_spectrum

    if (verbose >= 3) print*,'Start of ', nameOfRoutine

      ! first put default values
      call settings_fill_default()
  
      if (namelist_file == "None") then
        if (verbose >= 3) print*,'No namelist file to read!', namelist_file
        errorstatus = success
        return
      end if

  

      ! read name list parameter file
      open(7, file=namelist_file,delim='APOSTROPHE')
      read(7,nml=inoutput_mode)
      read(7,nml=output)
      read(7,nml=run_mode)
      read(7,nml=surface_params)
      read(7,nml=gas_abs_mod)
      read(7,nml=hyd_opts)
      read(7,nml=moments)
      read(7,nml=radar_simulator)

      close(7)

      call add_settings(err)
      if (err /= 0) then
          msg = 'error in add_settings!'
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
      end if

      call test_settings(err)
      if (err /= 0) then
          msg = 'error in test_settings!'
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
      end if

    errorstatus = err
    if (verbose >= 3) print*,'End of ', nameOfRoutine
    return
  end subroutine settings_read
    
  subroutine test_settings(errorstatus)
    use kinds
    implicit none
    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=15) :: nameOfRoutine = 'test_settings'
    !test for settings go here
    if (verbose >= 4) print*,'Start of ', nameOfRoutine
    err = 0
    call assert_false(err,MOD(radar_nfft, 2) == 1,&
        "radar_nfft has to be even") 
    call assert_true(err,(gas_mod == "L93") .or. (gas_mod == "R98"),&
        "gas_mod has to be L93 or R98") 
    if (hydro_fullSpec) then
      call assert_true(err,in_python,&
          "hydro_fullSpec works only in python!") 
    end if
    if (jacobian_mode) then
      call assert_true(err,randomseed/=0,&
          "randomniness not allowed in jacobian mode") 
    end if

    call assert_true(err,(radar_nPeaks == 1),&
        "radar_nPeaks higher than one not implemented yet!") 

    if (err /= 0) then
      msg = 'value in settings not allowed'
      call report(err, msg, nameOfRoutine)
      errorstatus = err
      return
    end if

    errorstatus = err
    if (verbose >= 4) print*,'End of ', nameOfRoutine
    return

  end subroutine test_settings

  subroutine add_settings(errorstatus)
    use kinds
    implicit none
    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=15) :: nameOfRoutine = 'add_settings'
    !additional variables derived from others go here

    if (verbose >= 4) print*,'Start of ', nameOfRoutine
    !mix some variables to make new ones:
    radar_nfft_aliased = radar_nfft *(1+2*radar_aliasing_nyquist_interv)
    radar_maxTurbTerms = radar_nfft_aliased * 12

    !in python some options are missing:
    if (in_python) add_obs_height_to_layer = .false.
 
    errorstatus = err
    if (verbose > 3) call print_settings()
    if (verbose >= 4) print*,'End of ', nameOfRoutine
    return

  end subroutine add_settings

  subroutine settings_fill_default
    use kinds
    implicit none

    character(len=14) :: nameOfRoutine = 'settings_fill_default'
    
    if (verbose >= 2) print*,'Start of ', nameOfRoutine

    
        hydro_threshold = 1.d-10   ! [kg/kg] 


        !set namelist defaults!
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
        save_psd=.false.
        ! sec output
        obs_height=833000.
        units='T'
        outpol='VH'
        freq_str=''
        file_desc=''
        creator='Pamtrauser'
        add_obs_height_to_layer = .false.
        ! sec run_mode
        active=.true.
        passive=.true.
        radar_mode="simple" !"splitted"|"moments"|"spectrum"
        randomseed = 0
        ! sec surface params
        ground_type='L'
        salinity=33.0
        emissivity=0.6
        ! sec gas_abs_mod
        lgas_extinction=.true.
        gas_mod='R98'
        ! sec hyd_opts
        lhyd_extinction=.true.
        lphase_flag = .true.
        hydro_fullSpec = .false.
        hydro_limit_density_area = .true.
        hydro_softsphere_min_density = 10. !kg/m^3
!        ! sec moments
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
        !radar noise at 1km in same unit as Ze 10*log10(mm⁶/m³). noise is calculated with noise = radar_pnoise0 + 20*log10(range/1000)
        radar_pnoise0=-32.23 ! mean value for BArrow MMCR during ISDAC

        radar_airmotion = .false.
        radar_airmotion_model = "step" !"constant","linear","step"
        radar_airmotion_vmin = -4.d0
        radar_airmotion_vmax = +4.d0
        radar_airmotion_linear_steps = 30
        radar_airmotion_step_vmin = 0.5d0

        radar_aliasing_nyquist_interv = 1
        radar_save_noise_corrected_spectra = .false.
        radar_use_hildebrand = .false.
        radar_min_spectral_snr = 1.2!threshold for peak detection. if radar_no_Ave >> 150, it can be set to 1.1
        radar_convolution_fft = .true. !use fft for convolution of spectrum. is alomst 10 times faster, but can introduce aretfacts for radars with *extremely* low noise levels or if noise is turned off at all.
        radar_K2 = 0.93 ! dielectric constant |K|² (always for liquid water by convention) for the radar equation
        radar_noise_distance_factor = 1.25
        radar_receiver_uncertainty_std = 0.d0 !dB
        radar_nPeaks = 1 !number of peaks the radar simulator is looking for
        radar_smooth_spectrum = .true.

        ! create frequency string if not set in pamtra
        if (freq_str == "") then
             ! get integer and character frequencies
            frq_str_s = "_"//frqs_str(1)
            if (nfrq .eq. 1) then
                frq_str_e = ""
            else
                frq_str_e = "-"//frqs_str(nfrq)
            end if
            freq_str = frq_str_s//frq_str_e
            
        end if 
    
    
    
    end subroutine settings_fill_default
    
    !for debuging
    subroutine print_settings()

      print*, 'add_obs_height_to_layer: ', add_obs_height_to_layer
      print*, 'jacobian_mode: ', jacobian_mode
      print*, 'radar_nfft: ', radar_nfft
      print*, 'input_path: ', input_path
      print*, 'creator: ', creator
      print*, 'radar_mode: ', radar_mode
      print*, 'lphase_flag: ', lphase_flag
      print*, 'radar_pnoise0: ', radar_pnoise0
      print*, 'data_path: ', data_path
      print*, 'radar_min_spectral_snr: ', radar_min_spectral_snr
      print*, 'radar_aliasing_nyquist_interv: ', radar_aliasing_nyquist_interv
      print*, 'radar_airmotion_linear_steps: ', radar_airmotion_linear_steps
      print*, 'radar_airmotion: ', radar_airmotion
      print*, 'ground_type: ', ground_type
      print*, 'obs_height: ', obs_height
      print*, 'crm_case: ', crm_case
      print*, 'write_nc: ', write_nc
      print*, 'outpol: ', outpol
      print*, 'radar_no_ave: ', radar_no_ave
      print*, 'input_type: ', input_type
      print*, 'dump_to_file: ', dump_to_file
      print*, 'passive: ', passive
      print*, 'radar_airmotion_model: ', radar_airmotion_model
      print*, 'crm_data: ', crm_data
      print*, 'tmp_path: ', tmp_path
      print*, 'lgas_extinction: ', lgas_extinction
      print*, 'crm_constants: ', crm_constants
      print*, 'units: ', units
      print*, 'gas_mod: ', gas_mod
      print*, 'moments_file: ', moments_file
      print*, 'radar_receiver_uncertainty_std: ', radar_receiver_uncertainty_std
      print*, 'hydro_threshold: ', hydro_threshold
      print*, 'hydro_fullSpec: ', hydro_fullSpec
      print*, 'hydro_limit_density_area: ', hydro_limit_density_area
      print*, 'hydro_softsphere_min_density: ', hydro_softsphere_min_density
      print*, 'radar_noise_distance_factor: ', radar_noise_distance_factor
      print*, 'radar_airmotion_step_vmin: ', radar_airmotion_step_vmin
      print*, 'crm_data2: ', crm_data2
      print*, 'radar_use_hildebrand: ', radar_use_hildebrand
      print*, 'radar_convolution_fft: ', radar_convolution_fft
      print*, 'radar_smooth_spectrum', radar_smooth_spectrum
      print*, 'active: ', active
      print*, 'radar_max_v: ', radar_max_v
      print*, 'radar_save_noise_corrected_spectra: ', radar_save_noise_corrected_spectra
      print*, 'n_moments: ', n_moments
      print*, 'lhyd_extinction: ', lhyd_extinction
      print*, 'radar_k2: ', radar_k2
      print*, 'radar_nPeaks', radar_nPeaks
      print*, 'salinity: ', salinity
      print*, 'radar_airmotion_vmax: ', radar_airmotion_vmax
      print*, 'output_path: ', output_path
      print*, 'radar_min_v: ', radar_min_v
      print*, 'freq_str: ', freq_str
      print*, 'file_desc: ', file_desc
      print*, 'radar_airmotion_vmin: ', radar_airmotion_vmin
      print*, 'emissivity: ', emissivity
      print*, 'save_psd: ', save_psd
      print*, "randomseed", randomseed

    end subroutine print_settings
    
end module settings

module nml_params
  ! Description:
  ! Definition of all name list paramters for pamtra
  !and global settings!


  use kinds

  implicit none

  !!Global Stettings
  integer, parameter :: MAXV = 64,   &
       MAXLAY = 200, &
       maxleg = 200, &
       maxfreq = 40, &
       nummu = 16, & ! no. of observation angles
       NSTOKES = 2, &
       NOUTLEVELS = 2
  integer, parameter :: SRC_CODE = 2,&
       NUMAZIMUTHS =1,&
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
  real(kind=dbl) :: snow_density, graupel_density, hail_density
  real(kind=dbl) :: salinity         ! sea surface salinity

  logical :: dump_to_file   ! flag for profile and ssp dump
  logical :: lphase_flag, &        ! flag for phase function calculation
       lgas_extinction, &    ! gas extinction desired
       lhyd_extinction, &    ! hydrometeor extinction desired
       write_nc, &	   ! write netcdf output
       active, &  	   ! calculate active stuff
       passive		   ! calculate passive stuff (with RT3)

  character(5) :: EM_ice, EM_snow, EM_grau, EM_hail
  character(1) :: SD_cloud, SD_ice, SD_rain, SD_snow, SD_grau, SD_hail
  character(3) :: gas_mod
  character(20) :: moments_file,file_desc
  character(100) :: input_path, output_path, tmp_path,creator, data_path
  character(13) :: freq_str
  character(2) :: OUTPOL
  character(1) :: GROUND_TYPE, UNITS



contains

  subroutine nml_params_read(namelist_file)

    character(300), intent(in) ::namelist_file

    ! name list declarations

    namelist / verbose_mode / verbose
    namelist / inoutput_mode / write_nc, input_path, output_path,&
         tmp_path, dump_to_file, data_path
    namelist / output / obs_height,units,outpol,freq_str,file_desc,creator
    namelist / run_mode / active, passive
    namelist / surface_params / ground_type,salinity, emissivity
    namelist / gas_abs_mod / lgas_extinction, gas_mod
    namelist / hyd_opts / lhyd_extinction, lphase_flag
    namelist / cloud_params / SD_cloud
    namelist / ice_params / SD_ice, EM_ice
    namelist / rain_params / SD_rain, N_0rainD
    namelist / snow_params / SD_snow, N_0snowDsnow, EM_snow, snow_density, SP, isnow_n0, liu_type
    namelist / graupel_params / SD_grau, N_0grauDgrau, EM_grau, graupel_density
    namelist / hail_params / SD_hail, N_0hailDhail, EM_hail, hail_density
    namelist / moments / n_moments, moments_file


    !set namelist defaults!
    verbose=0

    write_nc=.true.
    dump_to_file=.false.
    input_path='profile/'
    output_path='output/'
    tmp_path='/tmp/'
    data_path='data/'

    obs_height=833000.
    units='T'
    outpol='VH'
    freq_str=''
    file_desc=''
    creator='Pamtrauser'

    active=.true.
    passive=.true.

    ground_type='S'
    salinity=33.0
    emissivity=0.6

    lgas_extinction=.true.
    gas_mod='R98'

    lhyd_extinction=.true.
    lphase_flag = .true.

 	SD_cloud='C'

 	SD_ice='C'
    EM_ice='mieic'

    SD_rain='C'
    N_0rainD=8.0

    SD_snow='C'
    N_0snowDsnow=7.628 
    EM_snow='densi'
    snow_density=200.d0
    SP=0.2 
    isnow_n0=1
    liu_type=8

    SD_grau='C'
    N_0grauDgrau=4.0 
    EM_grau='densi'
    graupel_density=400.d0

    SD_hail='C'
    N_0hailDhail=4.0
    EM_hail='densi'
    hail_density=917.d0

    n_moments=1
    moments_file='snowCRYSTAL'

    ! read name list parameter file

    open(7, file=namelist_file,delim='APOSTROPHE')
    read(7,nml=verbose_mode)
    read(7,nml=inoutput_mode)
    read(7,nml=output)
    read(7,nml=run_mode)
    read(7,nml=surface_params)
    read(7,nml=gas_abs_mod)
    read(7,nml=hyd_opts)
    read(7,nml=cloud_params)
    read(7,nml=ice_params)
    read(7,nml=rain_params)
    read(7,nml=snow_params)
    read(7,nml=graupel_params)
    read(7,nml=hail_params)
    read(7,nml=moments)
    close(7)

    if (n_moments .ne. 1 .and. n_moments .ne. 2) stop "n_moments is not 1 or 2"

    if (verbose .gt. 1) print *,"PASSIVE: ", passive, "ACTIVE: ", active

    return
  end subroutine nml_params_read
end module nml_params

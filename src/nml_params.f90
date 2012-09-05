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
       zeSplitUp     ! save Ze and Att for every hydrometeor seperately. has only effect on netcdf file!

  character(5) :: EM_ice, EM_snow, EM_grau, EM_hail
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


contains

  subroutine nml_params_read

    use file_mod, only: namelist_file

    ! name list declarations
    namelist / verbose_mode / verbose
    namelist / inoutput_mode / input_path, output_path,&
         tmp_path, dump_to_file, write_nc, data_path,&
         input_type, crm_case, crm_data, crm_data2, crm_constants
    namelist / output / obs_height,units,outpol,freq_str,file_desc,creator,zeSplitUp
    namelist / run_mode / active, passive,rt_mode
    namelist / surface_params / ground_type,salinity, emissivity
    namelist / gas_abs_mod / lgas_extinction, gas_mod
    namelist / hyd_opts / lhyd_extinction, lphase_flag
    namelist / cloud_params / SD_cloud
    namelist / ice_params / SD_ice, EM_ice
    namelist / rain_params / SD_rain, N_0rainD, use_rain_db
    namelist / snow_params / SD_snow, N_0snowDsnow, EM_snow, use_snow_db, as_ratio,snow_density, SP, isnow_n0, liu_type
    namelist / graupel_params / SD_grau, N_0grauDgrau, EM_grau, graupel_density
    namelist / hail_params / SD_hail, N_0hailDhail, EM_hail, hail_density
    namelist / moments / n_moments, moments_file



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
    ! sec output
    obs_height=833000.
    units='T'
    outpol='VH'
    freq_str=''
    file_desc=''
    creator='Pamtrauser'
    zeSplitUp = .true.    
    ! sec run_mode
    active=.true.
    passive=.true.
    rt_mode='rt3'
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
    ! sec ice_params
    SD_ice='C'
    EM_ice='mieic'
    ! sec rain_params
    SD_rain='C'
    N_0rainD=8.0
    use_rain_db=.false.
    ! sec snow_params
    SD_snow='C'
    N_0snowDsnow=7.628 
    EM_snow='densi'
    use_snow_db=.false.
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
    close(7)

    if (n_moments .ne. 1 .and. n_moments .ne. 2) stop "n_moments is not 1 or 2"

    if (verbose .gt. 1) print *,"PASSIVE: ", passive, "ACTIVE: ", active


    return
  end subroutine nml_params_read
end module nml_params

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




  integer :: verbose, n_moments, isnow_n0, liu_type

  real(kind=dbl) :: obs_height     ! upper level output height [m] (> 100000. for satellite)

  real(kind=dbl) :: emissivity

  real(kind=dbl) :: N_0snowDsnow, N_0grauDgrau, N_0rainD, SP

  logical :: dump_to_file   ! flag for profile and ssp dump

  logical :: lphase_flag, &        ! flag for phase function calculation
	     lgas_extinction, &    ! gas extinction desired
	     lhyd_extinction, &    ! hydrometeor extinction desired
	     write_nc, &	   ! write netcdf output
	     active, &  	   ! calculate active stuff
	     passive		   ! calculate passive stuff (with RT3)

  character(5) :: EM_snow, EM_grau, EM_ice
  character(3) :: SD_snow, SD_grau, SD_rain

  character(3) :: gas_mod

  character(20) :: moments_file

  character(100) :: input_path, output_path, tmp_path,creator, data_path

  character :: OUTPOL*2, GROUND_TYPE*1, UNITS*1

  real(kind=dbl) :: salinity         ! sea surface salinity

  contains
    
    subroutine nml_params_read(namelist_file)

    character(300), intent(in) ::namelist_file

      ! name list declarations

      namelist / verbose_mode / verbose
      namelist / inoutput_mode / write_nc, input_path, output_path,&
         tmp_path, dump_to_file, data_path
      namelist / output / obs_height,units,outpol,creator
      namelist / run_mode / active, passive
      namelist / surface_params / ground_type,salinity, emissivity
      namelist / gas_abs_mod / lgas_extinction, gas_mod
      namelist / hyd_opts / lhyd_extinction, lphase_flag
      namelist / snow_params / SD_snow, N_0snowDsnow, EM_snow, SP, isnow_n0, liu_type
      namelist / graupel_params / SD_grau, N_0grauDgrau, EM_grau
      namelist / ice_params / EM_ice
      namelist / rain_params / SD_rain, N_0rainD
      namelist / moments / n_moments, moments_file


      !set namelist defaults!
      verbose=0

      write_nc=.true.
      dump_to_file=.false.
      input_path='profile'
      output_path='output'
      tmp_path='/tmp/'
      data_path='/home/mech/models/pamtra/data/'

      obs_height=833000.
      units='T'
      outpol='VH' 
      creator='Pamtra'

      active=.true.
      passive=.true.

      ground_type='S'
      salinity=33.0
      emissivity=0.6

      lgas_extinction=.true.
      gas_mod='R98'

      lhyd_extinction=.true.
      lphase_flag = .true.

      SD_snow='Exp' 
      N_0snowDsnow=7.628 
      EM_snow='icesf' 
      SP=0.2 
      isnow_n0=1
      liu_type=8

      SD_grau='Exp' 
      N_0grauDgrau=4.0 
      EM_grau='surus'

      EM_ice='mieic'

      SD_rain='Exp' 
      N_0rainD=8.0

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
      read(7,nml=snow_params)
      read(7,nml=graupel_params)
      read(7,nml=ice_params)
      read(7,nml=rain_params)
      read(7,nml=moments)
      close(7)

      if (n_moments .ne. 1 .and. n_moments .ne. 2) stop "n_moments is not 1 or 2"

      if (verbose .gt. 1) print *,"PASSIVE: ", passive, "ACTIVE: ", active

      return
    end subroutine nml_params_read
end module nml_params

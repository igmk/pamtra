subroutine pyPamtraLib(input_file,frequency&
, & !settings
set_verbose,set_write_nc,set_dump_to_file,set_input_path,set_output_path,set_tmp_path,&
set_data_path,set_obs_height,set_units,set_outpol,set_freq_str,set_file_desc,set_creator,&
set_active,set_passive,set_ground_type,set_salinity,set_emissivity,set_lgas_extinction,&
set_gas_mod,set_lhyd_extinction,set_lphase_flag,set_SD_snow,set_N_0snowDsnow,set_EM_snow,&
set_SP,set_isnow_n0,set_liu_type,set_SD_grau,set_N_0grauDgrau,set_EM_grau,set_EM_ice,set_SD_rain,&
set_N_0rainD,set_n_moments,set_moments_file&
,& !meta out
gitVersion,gitHash &
)

! ,&
! in_year,in_month,in_day,in_time,in_ngridx,in_ngridy,in_nlyr,in_deltax,in_deltay,in_&
! in_lat,in_lon,in_lfrac,in_relhum_lev,in_press_lev,in_temp_lev,in_hgt_lev,in_&
! in_model_i,in_model_j,in_wind10u,in_wind10v,in_iwv,in_cwp,in_&
! in_iwp,in_rwp,in_swp,in_gwp,in_hwp,in_cwc_q,in_iwc_q,in_rwc_q,in_swc_q,in_gwc_q)

  use kinds
  use constants !physical constants live here
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use vars_profile
  use double_moments_module !double moments variables are stored here
  use mod_io_strings !some strings for nice filenames




  !     The code reads a full (e.g. COSMO) grid and computes for each 
  !     profile the radiative transfer for the given frequencies                       
  !                                                                      
  !     By convention, the quantities followed  by "_lev"
  !     are given at the layer heights while the quantitites w/o
  !     "_lev" are layer average quantities                           
  !                                                                       

  implicit none

!!! internal "handle command line parameters" !!! 

  integer :: inarg, ff
  character(40),intent(out) :: gitHash, gitVersion
  character(6) :: formatted_frqstr !function call

!!! set by "handle command line parameters" !!! 

  character(99),intent(in)  :: input_file !name of profile
  real(kind=sgl), intent(in) :: frequency
  real(kind=sgl) :: freq
  !!Set by namelist file
  integer :: set_verbose, set_n_moments, set_isnow_n0, set_liu_type

  real(kind=sgl) :: set_obs_height     ! upper level output height [m] (> 100000. for satellite)
  real(kind=sgl) :: set_emissivity
  real(kind=sgl) :: set_N_0snowDsnow, set_N_0grauDgrau, set_N_0rainD, set_SP
  real(kind=sgl) :: set_salinity         ! sea surface salinity

  logical :: set_dump_to_file   ! flag for profile and ssp dump
  logical :: set_lphase_flag, &        ! flag for phase function calculation
       set_lgas_extinction, &    ! gas extinction desired
       set_lhyd_extinction, &    ! hydrometeor extinction desired
       set_write_nc, &	   ! write netcdf output
       set_active, &  	   ! calculate active stuff
       set_passive		   ! calculate passive stuff (with RT3)

  character(5) :: set_EM_snow, set_EM_grau, set_EM_ice
  character(3) :: set_SD_snow, set_SD_grau, set_SD_rain, set_gas_mod
  character(20) :: set_moments_file,set_file_desc
  character(100) :: set_input_path, set_output_path, set_tmp_path,set_creator, set_data_path
  character(13) :: set_freq_str
  character(2) :: set_OUTPOL
  character(1) :: set_GROUND_TYPE, set_UNITS




!f2py intent(in) :: input_file,frequency
!settings
!f2py intent(in) :: set_verbose,set_write_nc,set_dump_to_file,set_input_path,set_output_path,set_tmp_path
!f2py intent(in) :: set_data_path,set_obs_height,set_units,set_outpol,set_freq_str,set_file_desc,set_creator
!f2py intent(in) :: set_active,set_passive,set_ground_type,set_salinity,set_emissivity,set_lgas_extinction
!f2py intent(in) :: set_gas_mod,set_lhyd_extinction,set_lphase_flag,set_SD_snow,set_N_0snowDsnow,set_EM_snow
!f2py intent(in) :: set_SP,set_isnow_n0,set_liu_type,set_SD_grau,set_N_0grauDgrau,set_EM_grau,set_EM_ice,set_SD_rain
!f2py intent(in) :: set_N_0rainD,set_n_moments,set_moments_file
!meta out
!f2py intent(out) :: gitVersion,gitHash


!!!loop variables
  integer ::  fi,nx, ny

!!!output variables
  character(300) ::nc_out_file

print *,"Hello"
print *,input_file,frequency


  !get git data
  call versionNumber(gitVersion,gitHash)

print *,gitVersion,gitHash

nfrq = 1
fi = 1


print *,frequency

!load settings, uggly but neccessary!
verbose = set_verbose
write_nc = set_write_nc
dump_to_file = set_dump_to_file
input_path = set_input_path
output_path = set_output_path
tmp_path = set_tmp_path
data_path = set_data_path
obs_height = set_obs_height
units = set_units
outpol = set_outpol
freq_str = set_freq_str
file_desc = set_file_desc
creator = set_creator
active = set_active
passive = set_passive
ground_type = set_ground_type
salinity = set_salinity
emissivity = set_emissivity
lgas_extinction = set_lgas_extinction
gas_mod = set_gas_mod
lhyd_extinction = set_lhyd_extinction
lphase_flag = set_lphase_flag
SD_snow = set_SD_snow
N_0snowDsnow = set_N_0snowDsnow
EM_snow = set_EM_snow
SP = set_SP
isnow_n0 = set_isnow_n0
liu_type = set_liu_type
SD_grau = set_SD_grau
N_0grauDgrau = set_N_0grauDgrau
EM_grau = set_EM_grau
EM_ice = set_EM_ice
SD_rain = set_SD_rain
N_0rainD = set_N_0rainD
n_moments = set_n_moments
moments_file = set_moments_file



  freq = frequency

  ! create frequency string of not set in pamtra
  if (freq_str .eq. "") then
     read(freq_str,"(f3.2)") freq
  end if

print*, frequency, freq_str

  if (verbose .gt. 1) print *,"input_file: ",input_file(:len_trim(input_file)),&
       " freq: ",freq_str

!!! read n-moments file
  if (n_moments .eq. 2) call double_moments_module_read(moments_file) !from double_moments_module.f90

!!! read the data
  call vars_profile_read_profile(input_file) !from vars_atmosphere.f90

year = profiles_year
month = profiles_month
day = profiles_day
time = profiles_time
ngridx = profiles_ngridx
ngridy = profiles_ngridy
nlyr = profiles_nlyr
deltax = profiles_deltax
deltay = profiles_deltay
  ! now allocate variables
  call allocate_vars


  if (write_nc .eqv. .false.) call mod_io_strings_get_filename()


  if (verbose .gt. 1) print*, 'Start loop over frequencies & profiles!'



     grid_y: do ny = 1, ngridy !ny_in, ny_fin  
        grid_x: do nx = 1, ngridx !nx_in, nx_fin   

         !   ground_temp = profiles(nx,ny)%temp_lev(0)       ! K
         lat = profiles(nx,ny)%latitude                  ! °
         lon = profiles(nx,ny)%longitude                 ! °
         lfrac = profiles(nx,ny)%land_fraction
         relhum_lev = profiles(nx,ny)%relhum_lev         ! %
         press_lev = profiles(nx,ny)%press_lev           ! Pa
         temp_lev = profiles(nx,ny)%temp_lev             ! K
         hgt_lev = profiles(nx,ny)%hgt_lev               ! m

         model_i = profiles(nx,ny)%isamp
         model_j = profiles(nx,ny)%jsamp
         wind10u = profiles(nx,ny)%wind_10u
         wind10v = profiles(nx,ny)%wind_10v

         iwv = profiles(nx,ny)%iwv
         cwp = profiles(nx,ny)%cwp
         iwp = profiles(nx,ny)%iwp
         rwp = profiles(nx,ny)%rwp
         swp = profiles(nx,ny)%swp
         gwp = profiles(nx,ny)%gwp
         hwp = profiles(nx,ny)%hwp


         cwc_q = profiles(nx,ny)%cloud_water_q           ! kg/kg
         iwc_q = profiles(nx,ny)%cloud_ice_q             ! kg/kg
         rwc_q = profiles(nx,ny)%rain_q                  ! kg/kg
         swc_q = profiles(nx,ny)%snow_q                  ! kg/kg
         gwc_q = profiles(nx,ny)%graupel_q               ! kg/kg

         if (n_moments .eq. 2) then
            hwc_q = profiles(nx,ny)%hail_q              ! kg/kg
            cwc_n = profiles(nx,ny)%cloud_water_n       ! #/kg
            iwc_n = profiles(nx,ny)%cloud_ice_n         ! #/kg
            rwc_n = profiles(nx,ny)%rain_n              ! #/kg
            swc_n = profiles(nx,ny)%snow_n              ! #/kg
            gwc_n = profiles(nx,ny)%graupel_n           ! #/kg
            hwc_n = profiles(nx,ny)%hail_n              ! #/kg
         end if

           !run the model
           call run_rt3(nx,ny,fi,frequency,freq_str)

        end do grid_x
     end do grid_y

end subroutine pyPamtraLib
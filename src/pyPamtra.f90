subroutine pyPamtra(input_file,namelist_file,frequency,&
in_year,in_month,in_day,in_time,in_ngridx,in_ngridy,in_nlyr,in_deltax,in_deltay,in_&
in_lat,in_lon,in_lfrac,in_relhum_lev,in_press_lev,in_temp_lev,in_hgt_lev,in_&
in_model_i,in_model_j,in_wind10u,in_wind10v,in_iwv,in_cwp,in_&
in_iwp,in_rwp,in_swp,in_gwp,in_hwp,in_cwc_q,in_iwc_q,in_rwc_q,in_swc_q,in_gwc_q)






  use kinds
  use constants !physical constants live here
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
!   use vars_profile
  use double_moments_module !double moments variables are stored here
  use mod_io_strings !some strings for nice filenames

  implicit none


!!! internal "handle command line parameters" !!! 

  integer :: inarg, ff
  character(40) :: gitHash, gitVersion
  character(6) :: formatted_frqstr !function call

!!! set by "handle command line parameters" !!! 

  character(99)  :: input_file !name of profile
  character(300) :: namelist_file
  character(6), dimension(maxfreq) :: frqs_str !from commandline
  character(7) :: frq_str_s,frq_str_e

!!!loop variables
  integer ::  fi,nx, ny

!!!output variables
  character(300) ::nc_out_file



!!! read variables from namelist file
call nml_params_read(namelist_file) !from nml_params.f90


lat = in_latitude                  ! °
lon = in_longitude                 ! °
lfrac = in_land_fraction
relhum_lev = in_relhum_lev         ! %
press_lev = in_press_lev           ! Pa
temp_lev = in_temp_lev             ! K
hgt_lev = in_hgt_lev               ! m

model_i = in_isamp
model_j = in_jsamp
wind10u = in_wind_10u
wind10v = in_wind_10v

iwv = in_iwv
cwp = in_cwp
iwp = in_iwp
rwp = in_rwp
swp = in_swp
gwp = in_gwp
hwp = in_hwp


cwc_q = in_cloud_water_q           ! kg/kg
iwc_q = in_cloud_ice_q             ! kg/kg
rwc_q = in_rain_q                  ! kg/kg
swc_q = in_snow_q                  ! kg/kg
gwc_q = in_graupel_q               ! kg/kg



nx = 1
ny = 1
fi = 1




!run the model
call run_rt3(nx,ny,fi,freqs(fi),frqs_str(fi))

end subroutine pyPamtra
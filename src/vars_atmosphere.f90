module vars_atmosphere

  use kinds
  use settings, only: maxfreq

  implicit none
  save

  integer(kind=long) :: nlyr
  integer(kind=long) :: ngridx, ngridy
  character(2) :: month, day
  character(4) :: year, time
  character(12) :: date_str
  real(kind=sgl) :: deltax, deltay

  integer(kind=long), allocatable, dimension(:) :: rt_nlegen

  !is allocated in pamtra.f90!

  real(kind=dbl), allocatable, dimension(:) :: relhum_lev,&
       press_lev, &
       temp_lev, &
       hgt_lev

  real(kind=dbl), allocatable, dimension(:) :: press, &
       temp,&
       relhum,&
       vapor_pressure, &
       rho_vap, &
       q_hum,&
       hgt

  real(kind=dbl), allocatable, dimension(:) :: cwc_q, &
       iwc_q, &
       rwc_q, &
       swc_q, &
       gwc_q, &
       hwc_q

  real(kind=dbl), allocatable, dimension(:,:) :: q_hydro


  real(kind=dbl), allocatable, dimension(:) :: cwc_n, &
       iwc_n, &
       rwc_n, &
       swc_n, &
       gwc_n, &
       hwc_n

  real(kind=dbl), allocatable, dimension(:) :: delta_hgt_lev


  integer(kind=long) :: alloc_status

  integer(kind=long) :: model_i, model_j
  real(kind=sgl) :: lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp,hwp

end module vars_atmosphere

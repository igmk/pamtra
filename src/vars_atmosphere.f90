module vars_atmosphere

  use kinds
  use nml_params, only: maxfreq
  implicit none
  save

  integer :: nlyr, nfrq
  integer :: ngridx, ngridy
  character(2) :: month, day
  character(4) :: year, time
  real(kind=sgl) :: deltax, deltay

  integer, allocatable, dimension(:) :: nlegen
  integer, allocatable, dimension(:) :: rt3nlegen

  !is allocated in pamtra.f90!
  real(kind=dbl), dimension(maxfreq) :: freqs

  real(kind=dbl), allocatable, dimension(:) :: relhum_lev,&
       press_lev, &
       temp_lev, &
       hgt_lev

  real(kind=dbl), allocatable, dimension(:) :: press, &
       temp,&
       relhum,&
       vapor_pressure, &
       rho_vap, &
       q_hum

  real(kind=dbl), allocatable, dimension(:) :: cwc_q, &
       iwc_q, &
       rwc_q, &
       swc_q, &
       gwc_q, &
       hwc_q

  real(kind=dbl), allocatable, dimension(:) ::	cwc_n, &
       iwc_n, &
       rwc_n, &
       swc_n, &
       gwc_n, &
       hwc_n

  real(kind=dbl), allocatable, dimension(:) :: kextatmo, &
       kexttot, &
       rt3kexttot,&
       salbtot, &
       rt3salbtot,&
       back, &
       g_coeff

  real(kind=dbl), allocatable, dimension(:,:) :: legen, &
       legen2, &
       legen3, &
       legen4

  real(kind=dbl), allocatable, dimension(:,:) :: rt3legen, &
       rt3legen2, &
       rt3legen3, &
       rt3legen4

  integer :: alloc_status

  integer :: model_i, model_j
  real(kind=sgl) :: lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp,hwp


  integer, allocatable, dimension(:,:) :: ics



end module vars_atmosphere

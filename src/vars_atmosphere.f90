module vars_atmosphere

  use kinds
  implicit none
  save

  type profile
     integer :: isamp
     integer :: jsamp
     real(kind=sgl) :: latitude
     real(kind=sgl) :: longitude
     real(kind=sgl) :: land_fraction
     real(kind=sgl) :: wind_10u, wind_10v
     real(kind=sgl) :: iwv
     real(kind=sgl) :: cwp
     real(kind=sgl) :: iwp
     real(kind=sgl) :: rwp
     real(kind=sgl) :: swp
     real(kind=sgl) :: gwp
     real(kind=sgl) :: hwp
     !  real(kind=dbl) :: ground_albedo
     !  real(kind=dbl) :: ground_index
     !  real(kind=dbl) :: ground_temperature
     real(kind=dbl), pointer :: hgt_lev(:)
     real(kind=dbl), pointer :: press_lev(:)
     real(kind=dbl), pointer :: temp_lev(:)
     real(kind=dbl), pointer :: relhum_lev(:)

     real(kind=dbl), pointer :: cloud_water_q(:)
     real(kind=dbl), pointer :: cloud_ice_q(:)
     real(kind=dbl), pointer :: rain_q(:)
     real(kind=dbl), pointer :: snow_q(:)
     real(kind=dbl), pointer :: graupel_q(:)
     real(kind=dbl), pointer :: hail_q(:)

     real(kind=dbl), pointer :: cloud_water_n(:)
     real(kind=dbl), pointer :: cloud_ice_n(:)
     real(kind=dbl), pointer :: rain_n(:)
     real(kind=dbl), pointer :: snow_n(:)
     real(kind=dbl), pointer :: graupel_n(:)
     real(kind=dbl), pointer :: hail_n(:)

     real(kind=dbl), pointer :: press(:)
     real(kind=dbl), pointer :: temp(:)
     real(kind=dbl), pointer :: relhum(:)
     real(kind=dbl), pointer :: vapor_pressure(:)
     real(kind=dbl), pointer :: rho_vap(:)
     real(kind=dbl), pointer :: q_hum(:)
  end type profile



  integer :: nlyr, nfrq
  integer :: ngridx, ngridy

  integer, allocatable, dimension(:) :: nlegen
  integer, allocatable, dimension(:) :: rt3nlegen

  !is allocated in pamtra.f90!
  real(kind=dbl), allocatable, dimension(:) :: freqs

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

  type(profile), allocatable :: profiles(:,:)
  integer, allocatable, dimension(:,:) :: ics

  character(2) :: month, day
  character(4) :: year, time
  real(kind=sgl) :: deltax, deltay

end module vars_atmosphere

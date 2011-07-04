module profile_type

use kinds

implicit none

type profile
  integer :: nlyr
  integer :: isamp
  integer :: jsamp
  integer, pointer :: nlegen(:)
  real :: latitude
  real :: longitude
  real :: land_fraction
  real :: wind_10u, wind_10v
  real :: iwv
  real :: cwp
  real :: iwp
  real :: rwp
  real :: swp
  real :: gwp
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

  real(kind=dbl), pointer :: press(:)
  real(kind=dbl), pointer :: temp(:)
  real(kind=dbl), pointer :: relhum(:)
  real(kind=dbl), pointer :: vapor_pressure(:)
  real(kind=dbl), pointer :: rho_vap(:)
  real(kind=dbl), pointer :: q_hum(:)

  real(kind=dbl), pointer :: kextatmo(:)
  real(kind=dbl), pointer :: kexttot(:)
  real(kind=dbl), pointer :: salbtot(:)
  real(kind=dbl), pointer :: g_coeff(:)
  real(kind=dbl), pointer :: back(:)
  real(kind=dbl), pointer :: legen(:,:)
  real(kind=dbl), pointer :: legen2(:,:)
  real(kind=dbl), pointer :: legen3(:,:)
  real(kind=dbl), pointer :: legen4(:,:)
end type profile

end module profile_type

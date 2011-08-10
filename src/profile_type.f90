module profile_type

use kinds

implicit none

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

end module profile_type

function spec2abs(spec_var,T,p,q)
! Converts specific variables (e.g. snow content in kg-snow per kg air)
! into absolute values (densities) (e.g. snow content in kg-snow per cubic meter air)
! by using the virtual temperature

  use kinds

  implicit none

  real(kind=dbl), parameter :: r_l = 287.05 ! gas constant of dry air

  real(kind=dbl), intent(in) :: spec_var,&! specific variable [kg/kg]
				T,&       ! air temperature [K]
				p,&       ! atmospheric pressure [Pa]
				q         ! specific humidity [kg/kg]


  real(kind=dbl) :: spec2abs ! variable in absolute values [kg/m^3]

  real(kind=dbl) :: T_v,&   ! virtual temperature [k]
		    rho_air ! density of moist air [kg/m^3]

  !calculate virtual temperature in K with air_temp and spec.humidity
  T_v = T + T * 0.6078 * q

  ! calculate density of moist air
  rho_air = p / (T_v * r_l)

  ! multiply [kg/kg] * [kg/m^3] = kg/m^3
  spec2abs = spec_var * rho_air

  return
end function spec2abs
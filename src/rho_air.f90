real(kind=dbl) function rho_air(T, P)

! This function returns the density of dry air in kg/m**3.
! The input parameters T and P are the temperature (K) and 
! pressure (pa)                           

!  use nml_params, only: verbose
  use kinds
        use report_module
  implicit none

  real(kind=dbl), intent (in) :: T, P
  real(kind=dbl) :: r = 28704.d-2

  rho_air = P / (T*r) 

  return 
end function rho_air

subroutine viscosity_air(T,eta)

! This function returns the dynamic viscosity of dry air in Pa s
! Sutherland law           
! coefficients from F. M. White, Viscous Fluid Flow, 2nd ed., McGraw-Hill, (1991). Kim et al., arXiv:physics/0410237v1
          

!  use settings, only: verbose
  use kinds
        use report_module
  implicit none

  real(kind=dbl), intent (in) :: T
  real(kind=dbl), intent (out) :: eta
  real(kind=dbl) :: mu0 = 1.716d-5 !Pas
  real(kind=dbl) :: T0 = 273.d0
  real(kind=dbl) :: C = 111.d0 !K
  eta = mu0 * ((T0+C)/(T+C)) * (T/T0)**1.5d0

  return 
end subroutine viscosity_air

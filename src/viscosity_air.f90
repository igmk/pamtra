real(kind=dbl) function viscosity_air(T)

! This function returns the dynamic viscosity of dry air in Pa s
           
! coefficients from F. M. White, Viscous Fluid Flow, 2nd ed., McGraw-Hill, (1991). Kim et al., arXiv:physics/0410237v1
          

!  use nml_params, only: verbose
  use kinds
        use report_module
  implicit none

  real(kind=dbl), intent (in) :: T
  real(kind=dbl) :: eta
  real(kind=dbl) :: mu0 = 18.27 * 1d-6 !Pas
  real(kind=dbl) :: T0 = 291.15 !K
  real(kind=dbl) :: C = 120 !K

  eta = mu0 * ((T0+C)/(T+C)) * (T/T0)**1.5d0

  viscosity_air = eta


  return 
end function viscosity_air

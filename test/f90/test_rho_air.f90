! Unit test for the stateless rho_air() ideal-gas-law function.
program test_rho_air
  use kinds, only: dbl

  implicit none

  real(kind=dbl) :: rho_air
  real(kind=dbl) :: got
  real(kind=dbl), parameter :: expected = 1.2923286909749199_dbl
  real(kind=dbl), parameter :: tol = 1.0d-9

  got = rho_air(273.15_dbl, 101325.0_dbl)

  if (abs(got - expected) > tol) then
    print *, "rho_air(273.15K, 101325Pa) mismatch. got=", got, " expected=", expected
    stop 1
  end if

  print *, "rho_air: OK", got
end program test_rho_air

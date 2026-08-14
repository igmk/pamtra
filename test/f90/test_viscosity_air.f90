! Unit test for the stateless viscosity_air() Sutherland-law subroutine.
program test_viscosity_air
  use kinds, only: dbl

  implicit none

  real(kind=dbl) :: got
  real(kind=dbl), parameter :: expected = 1.7167438768896050d-5
  real(kind=dbl), parameter :: tol = 1.0d-12

  call viscosity_air(273.15_dbl, got)

  if (abs(got - expected) > tol) then
    print *, "viscosity_air(273.15K) mismatch. got=", got, " expected=", expected
    stop 1
  end if

  print *, "viscosity_air: OK", got
end program test_viscosity_air
